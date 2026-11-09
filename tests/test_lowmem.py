"""Equivalence tests for debreak_lowmem (the --lowmem code path).

Every function in debreak_lowmem replaces an original one and must produce
byte-identical output. These tests run the original and the replacement on
the same inputs and compare; the end-to-end check of the whole pipeline in
both modes is in test_e2e.py (BOTH_MODES).
"""
import filecmp
import os
import random
import shutil
import subprocess
import sys

import pysam
import pytest

import debreak_detect as dd
import debreak_lowmem as lm
import debreak_merge_cluster as mc
import debreak_resdup_selfalignment as rd
from conftest import REPO

TEST_BAM = os.path.join(REPO, "testdata", "test_read.bam")
TRA_BAM = os.path.join(REPO, "testdata", "tra", "tra.bam")
TEST_REF = os.path.join(REPO, "testdata", "test_ref.fa")
LAUNCHER = os.path.join(REPO, "tests", "run_debreak.py")


def _tree(root):
    out = {}
    for d, _, files in os.walk(root):
        for f in files:
            p = os.path.join(d, f)
            out[os.path.relpath(p, root)] = p
    return out


def _assert_same_tree(a, b):
    ta, tb = _tree(a), _tree(b)
    assert sorted(ta) == sorted(tb)
    for rel in ta:
        assert filecmp.cmp(ta[rel], tb[rel], shallow=False), rel + " differs"
    return ta


# ---------------------------------------------------------------- slicing helpers

@pytest.mark.parametrize("seed", range(3))
def test_slice_bounds_matches_python_slicing(seed):
    rnd = random.Random(seed)
    for _ in range(3000):
        length = rnd.randint(0, 50)
        lo, hi = rnd.randint(-80, 80), rnd.randint(-80, 80)
        data = list(range(length))
        b = lm._slice_bounds(lo, hi, length)
        got = [] if b is None else data[b[0]:b[1]]
        assert got == data[lo:hi], (length, lo, hi)
        if b is not None:
            assert 0 <= b[0] < b[1] <= length


# ---------------------------------------------------------------- region detection

def _dense_regions(intervals, contiglength, maxdepth):
    """Region detection exactly as debreak_merge_cluster.cluster does it, with
    the per-base list. Reference implementation for lm._find_regions."""
    genomeposition = [0] * contiglength
    for lo, hi in intervals:
        original = genomeposition[lo:hi]
        genomeposition[lo:hi] = [mm + 1 for mm in original]
    svregion = []
    inblock = False
    threshold = 3
    for i in range(len(genomeposition)):
        if inblock:
            if genomeposition[i] >= max(maxdep / 10.0, threshold):
                localdep += [genomeposition[i]]
                if genomeposition[i] > maxdep:
                    maxdep = genomeposition[i]
            else:
                inblock = False
                end = i
                if maxdep <= maxdepth:
                    peakpos = localdep.index(maxdep)
                    peakleftsize = 0
                    for i in range(peakpos):
                        if localdep[peakpos - i - 1] >= maxdep / 10.0:
                            peakleftsize += 1
                        else:
                            break
                    svregion += [(start + peakpos - peakleftsize, end, maxdep)]
                start = 0
                end = 0
                maxdep = 0
        else:
            if genomeposition[i] > threshold:
                inblock = True
                localdep = [genomeposition[i]]
                start = i
                maxdep = genomeposition[i]
    return [c for c in svregion if c[2] < maxdepth]


def _random_intervals(rnd, contiglength, n):
    """Half-open 0-based intervals as cluster()/cluster_ins() derive them,
    including ones starting before the contig or running past its end."""
    out = []
    for _ in range(n):
        if rnd.random() < 0.2:  # a pile-up of overlapping intervals
            pos = rnd.randint(-120, contiglength + 50)
            size = rnd.randint(48, 400)
            for _ in range(rnd.randint(2, 60)):
                s = pos + rnd.randint(-30, 30)
                out.append((s - 1, s + size + rnd.randint(-5, 5) - 1))
        else:
            s = rnd.randint(-120, contiglength + 50)
            out.append((s - 1, s + rnd.randint(48, 3000) - 1))
    return out


@pytest.mark.parametrize("seed", range(6))
@pytest.mark.parametrize("contiglength", [150, 5_000, 60_000])
def test_find_regions_matches_dense_sweep(seed, contiglength):
    rnd = random.Random(seed)
    intervals = _random_intervals(rnd, contiglength, 200)
    for maxdepth in (10, 50):
        assert lm._find_regions(intervals, contiglength, maxdepth) == \
            _dense_regions(intervals, contiglength, maxdepth)


def test_find_regions_block_touching_contig_end_is_dropped_like_original():
    L = 1000
    at_end = [(900, L)] * 5  # depth 5 up to the last base: the original never closes it
    assert _dense_regions(at_end, L, 50) == []
    assert lm._find_regions(at_end, L, 50) == []
    before_end = [(900, L - 1)] * 5
    assert _dense_regions(before_end, L, 50) == [(900, 999, 5)]
    assert lm._find_regions(before_end, L, 50) == [(900, 999, 5)]


def test_find_regions_empty():
    assert lm._find_regions([], 1000, 50) == []
    assert lm._find_regions([(10, 20)] * 2, 1000, 50) == []  # depth 2 < threshold


# ---------------------------------------------------------------- cluster / cluster_ins

def _raw_lines(seed, contiglength, kind):
    """Raw SV signals as the driver hands them to cluster()/cluster_ins()."""
    rnd = random.Random(seed)
    lines = []
    rid = [0]

    def line(pos, size):
        rid[0] += 1
        return "\t".join(["chrT", str(pos), str(size), kind, f"r{rid[0]}_cigar1",
                          str(rnd.choice([0, 16])), str(rnd.choice([0, 0, 5, 20, 60, 60, 60]))])

    for _ in range(60):  # clusters with sizes on both sides of the 2000/3000 thresholds
        pos = rnd.randint(1, contiglength)
        size = rnd.choice([50, 100, 150, 300, 500, 1200, 1900, 2100, 2500, 3500, 6000])
        for _ in range(rnd.randint(1, 15)):
            lines.append(line(max(0, pos + rnd.randint(-40, 40)),
                              max(48, int(size * rnd.uniform(0.85, 1.15)))))
    for _ in range(3):  # pile-ups deeper than maxdepth
        pos = rnd.randint(1, contiglength)
        size = rnd.randint(60, 400)
        for _ in range(rnd.randint(60, 120)):
            lines.append(line(pos + rnd.randint(-10, 10), size + rnd.randint(-5, 5)))
    for _ in range(300):  # noise
        lines.append(line(rnd.randint(0, contiglength), rnd.randint(48, 3000)))
    for p in (0, 1, 50, 101, contiglength - 1, contiglength - 200, contiglength - 3000):
        for _ in range(5):  # contig edges
            lines.append(line(max(0, p + rnd.randint(0, 3)), rnd.choice([48, 250, 2990])))
    rnd.shuffle(lines)
    return lines


PARAMS = [(2, 20, None), (3, 50, None), (2, 30, 20), (4, 15, None)]


def _compare_cluster(tmp_path, lines, contiglength, mins, maxdepth, minmapq, ins):
    old, new = tmp_path / "old", tmp_path / "new"
    old.mkdir(parents=True)
    new.mkdir(parents=True)
    for outdir, mod in ((old, mc), (new, lm)):
        if ins:
            mod.cluster_ins(str(outdir) + "/", list(lines), "chrT", contiglength, mins, maxdepth, minmapq)
        else:
            for svtype in ("del", "dup", "inv"):
                mod.cluster(str(outdir) + "/", list(lines), "chrT", contiglength, mins, maxdepth,
                            svtype, minmapq)
    files = _assert_same_tree(old, new)
    return sum(1 for p in files.values() for _ in open(p))


@pytest.mark.parametrize("seed", range(3))
@pytest.mark.parametrize("mins, maxdepth, minmapq", PARAMS)
def test_cluster_matches_original(tmp_path, seed, mins, maxdepth, minmapq):
    lines = _raw_lines(seed, 100_000, "D-cigar")
    assert _compare_cluster(tmp_path, lines, 100_000, mins, maxdepth, minmapq, ins=False) > 0


@pytest.mark.parametrize("seed", range(3))
@pytest.mark.parametrize("mins, maxdepth, minmapq", PARAMS)
def test_cluster_ins_matches_original(tmp_path, seed, mins, maxdepth, minmapq):
    lines = _raw_lines(seed + 100, 100_000, "I-cigar")
    assert _compare_cluster(tmp_path, lines, 100_000, mins, maxdepth, minmapq, ins=True) > 0


@pytest.mark.parametrize("contiglength", [150, 5_000])
def test_cluster_matches_original_on_tiny_contigs(tmp_path, contiglength):
    # intervals before position 1 and past the contig end take the slicing edge cases
    _compare_cluster(tmp_path / "d", _raw_lines(7, contiglength, "D-cigar"), contiglength, 2, 20, None, ins=False)
    _compare_cluster(tmp_path / "i", _raw_lines(8, contiglength, "I-cigar"), contiglength, 2, 20, None, ins=True)


def test_cluster_empty_input(tmp_path):
    assert lm.cluster(str(tmp_path) + "/", [], "chr1", 1000, 3, 50, "del", None) == 0
    assert lm.cluster_ins(str(tmp_path) + "/", [], "chr1", 1000, 3, 50, None) == 0
    assert os.listdir(tmp_path) == []


# ---------------------------------------------------------------- split_raw_calls

def _driver_dicts(writepath, vcflist, chromosomes):
    """The raw-call dicts exactly as the driver builds them without --lowmem."""
    raw = {t: {c: [] for c in chromosomes} for t in ("ins", "del", "dup", "inv", "tra")}
    for filename in vcflist:
        try:
            allvcf = open(writepath + 'sv_raw_calls/' + filename, 'r').read().split('\n')[:-1]
            allvcf[0].split('\t')[0]
        except Exception:
            continue
        for event in allvcf:
            event_chrom = event.split('\t')[0]
            if 'I-' in event:
                event = event.split('\t')
                raw["ins"][event_chrom] += ['\t'.join(event[:7])]
            if 'D-' in event:
                raw["del"][event_chrom] += [event]
            if 'DUP-' in event:
                raw["dup"][event_chrom] += [event]
            if 'INV-' in event:
                raw["inv"][event_chrom] += [event]
            if 'TRA-' in event:
                e = event.split('\t')
                raw["tra"][event_chrom] += [e[0] + '\t' + e[1] + '\t' + e[2] + '\t' + e[3] + '\t' + e[7] + '\t' + e[5]]
    return raw


@pytest.mark.parametrize("golden", ["default", "tumor", "tra_default"])
def test_split_raw_calls_matches_driver_dicts(tmp_path, golden):
    lines = open(os.path.join(REPO, "tests", "golden", golden, "sv_raw_calls.txt")).read().splitlines()
    raw_dir = tmp_path / "sv_raw_calls"
    raw_dir.mkdir()
    # spread over several files like --samlist input, plus an empty and a missing file
    files = ["a.debreak.temp", "b.debreak.temp", "c.debreak.temp"]
    for i, f in enumerate(files):
        (raw_dir / f).write_text("".join(l + "\n" for l in lines[i::3]))
    (raw_dir / "empty.debreak.temp").write_text("")
    vcflist = files + ["empty.debreak.temp", "missing.debreak.temp"]
    chromosomes = sorted({l.split("\t")[0] for l in lines})
    wp = str(tmp_path) + "/"

    expected = _driver_dicts(wp, vcflist, chromosomes)
    splitdir = lm.split_raw_calls(wp, vcflist, chromosomes)
    assert splitdir == wp + "sv_raw_calls/lowmem_split/"
    total = 0
    for svtype in expected:
        for chrom in chromosomes:
            assert lm._read_split(splitdir, chrom, svtype) == expected[svtype][chrom], (svtype, chrom)
            total += len(expected[svtype][chrom])
    assert total == len(lines)
    lm.remove_split(splitdir)
    assert not os.path.exists(splitdir)
    assert sorted(os.listdir(raw_dir)) == sorted(files + ["empty.debreak.temp"])


def test_split_raw_calls_drops_insertion_sequences(tmp_path):
    (tmp_path / "sv_raw_calls").mkdir()
    (tmp_path / "sv_raw_calls" / "x.temp").write_text(
        "chr1\t100\t60\tI-cigar\tr1_cigar1\t0\t60\t500\tACGT\tTTTT\n"
        "chr1\t200\t70\tD-cigar\tr2_cigar1\t16\t55\n")
    splitdir = lm.split_raw_calls(str(tmp_path) + "/", ["x.temp"], ["chr1"])
    assert lm._read_split(splitdir, "chr1", "ins") == ["chr1\t100\t60\tI-cigar\tr1_cigar1\t0\t60"]
    assert lm._read_split(splitdir, "chr1", "del") == ["chr1\t200\t70\tD-cigar\tr2_cigar1\t16\t55"]
    assert lm._read_split(splitdir, "chr1", "tra") == []


# ---------------------------------------------------------------- detect_sortbam

def _detect(mod, bam, outdir, record_clip, rescuedup, tumor):
    for d in ("sv_raw_calls", "debreak_ins_workspace", "map_depth"):
        os.makedirs(os.path.join(outdir, d), exist_ok=True)
    chroms = list(pysam.AlignmentFile(bam).references)
    for chrom in chroms:
        mod.detect_sortbam(bam, outdir + "/", 48, 400_000_000, chrom, chroms,
                           record_clip, rescuedup, tumor)


@pytest.mark.parametrize("record_clip", [False, True])
@pytest.mark.parametrize("rescuedup", [False, True])
@pytest.mark.parametrize("tumor", [False, True])
def test_detect_sortbam_matches_original(tmp_path, record_clip, rescuedup, tumor):
    for name, bam in (("test_read", TEST_BAM), ("tra", TRA_BAM)):
        old, new = str(tmp_path / name / "old"), str(tmp_path / name / "new")
        _detect(dd, bam, old, record_clip, rescuedup, tumor)
        _detect(lm, bam, new, record_clip, rescuedup, tumor)
        files = _assert_same_tree(old, new)
        raw = [l for rel, p in files.items() if rel.startswith("sv_raw_calls") for l in open(p)]
        assert raw
        if name == "test_read" and not tumor:
            # split-read insertions carry the inserted sequence: that is what LazySeq fetches
            seg_ins = [l for l in raw if "I-segment" in l]
            assert seg_ins and all(len(l.rstrip("\n").split("\t")) == 8 and l.rstrip("\n").split("\t")[7]
                                   for l in seg_ins)
        if name == "tra":
            assert any("TRA-segment" in l for l in raw)
        if record_clip:
            assert any(rel.startswith("debreak_ins_workspace/readinfo_start_end_") for rel in files)


def test_lazyseq_fetches_the_primary_sequence_on_demand():
    with pysam.AlignmentFile(TEST_BAM) as bam, pysam.AlignmentFile(TEST_BAM) as seqfile:
        primary = next(a for a in bam.fetch("chr1")
                       if not a.is_secondary and not a.is_supplementary and a.has_tag("SA"))
        lazy = lm.LazySeq(seqfile, "chr1", primary.reference_start, primary.query_name)
        assert lazy.seq is None  # nothing read yet
        assert lazy[10:20] == primary.query_sequence[10:20]
        assert len(lazy) == len(primary.query_sequence)
        assert lazy[:] == primary.query_sequence


def test_lazyseq_unknown_read_is_empty():
    with pysam.AlignmentFile(TEST_BAM) as seqfile:
        lazy = lm.LazySeq(seqfile, "chr1", 0, "no_such_read")
        assert lazy[0:10] == ""
        assert len(lazy) == 0


# ---------------------------------------------------------------- DUP rescue

@pytest.fixture(scope="module")
def ref_copy(tmp_path_factory):
    """test_ref.fa copied out of the repo so that the .fai index lands in tmp."""
    dst = tmp_path_factory.mktemp("ref") / "test_ref.fa"
    shutil.copy(TEST_REF, dst)
    return str(dst)


def test_ref_slice_matches_string_slicing(ref_copy):
    seqs = {}
    for block in open(TEST_REF).read().split(">")[1:]:
        seqs[block.split("\n")[0].split(" ")[0]] = "".join(block.split("\n")[1:])
    assert seqs
    fasta = pysam.FastaFile(ref_copy)
    rnd = random.Random(0)
    for chrom, seq in seqs.items():
        L = len(seq)
        bounds = [(-5000, 100), (-10, 20), (0, 0), (0, 1), (L - 10, L + 500), (L, L + 1),
                  (L + 5, L + 10), (100, 50), (-L - 7, -L + 3)]
        bounds += [(rnd.randint(-L, 2 * L), rnd.randint(-L, 2 * L)) for _ in range(300)]
        for a, b in bounds:
            assert lm._ref_slice(fasta, chrom, a, b) == seq[a:b], (chrom, a, b)
    fasta.close()


@pytest.fixture(scope="module")
def default_run(tmp_path_factory):
    """Output directory of a default pipeline run (merged calls + raw calls)."""
    out = tmp_path_factory.mktemp("run") / "out"
    subprocess.run([sys.executable, LAUNCHER, "--bam", "testdata/test_read.bam",
                    "-o", str(out) + "/", "-t", "4"],
                   cwd=REPO, check=True, capture_output=True, stdin=subprocess.DEVNULL)
    return out


def test_identify_duplication_matches_original(default_run, ref_copy, tmp_path):
    # Without minimap2 on PATH the mapping step yields empty PAF files in both
    # variants; the reference-window extraction (the part that changed) is
    # exercised either way, and with minimap2 the rescue itself is compared too.
    vcflist = sorted(os.listdir(default_run / "sv_raw_calls"))
    outs = {}
    for tag, mod in (("old", rd), ("new", lm)):
        wp = tmp_path / tag
        shutil.copytree(default_run, wp)
        mod.identify_duplication(vcflist, str(wp) + "/", ref_copy)
        outs[tag] = wp
    files = _assert_same_tree(outs["old"], outs["new"])
    refseq = [rel for rel in files if rel.startswith("debreak_resdup_refseq/")]
    assert refseq  # one local reference window per merged insertion
    for rel in refseq:
        header, seq = open(files[rel]).read().split("\n")[:2]
        assert header.startswith(">chr1__") and len(seq) > 4000
