"""Characterization (golden master) tests for the whole `debreak` pipeline.

They pin the CURRENT behaviour, not necessarily the correct one, so that a
refactor can be checked for "no output changed". If a change is intended,
regenerate the snapshots with:  pytest --update-golden
"""
import os
import shutil
import subprocess
import sys

import pytest

from conftest import REPO

LAUNCHER = os.path.join(REPO, "tests", "run_debreak.py")
GOLDEN = os.path.join(REPO, "tests", "golden")
BAM = "testdata/test_read.bam"  # relative on purpose: it ends up in the VCF header
REF = "testdata/test_ref.fa"

# Header lines that legitimately differ between runs.
VOLATILE_HEADERS = ("##fileDate=", "##CommandLine=")


def _normalize_vcf(text):
    return "".join(
        line + "\n"
        for line in text.splitlines()
        if not line.startswith(VOLATILE_HEADERS)
    )


def _sorted_lines(text):
    return "".join(line + "\n" for line in sorted(text.splitlines()))


# Every pipeline scenario runs twice: with the original code path and with
# --lowmem (debreak_lowmem.py). Both must produce the same outputs, so they are
# checked against the same golden files. With --update-golden only the default
# run regenerates them; the --lowmem run is then compared against the fresh copy.
BOTH_MODES = pytest.mark.parametrize("lowmem", [False, True], ids=["default", "lowmem"])


def _run(outdir, *extra, bam=BAM, cwd=REPO, expect_rc=0, lowmem=False):
    """Run debreak. `bam=None` omits --bam (e.g. for --samlist runs)."""
    cmd = [sys.executable, LAUNCHER]
    if bam is not None:
        cmd += ["--bam", bam]
    if outdir is not None:
        cmd += ["-o", str(outdir) + "/"]
    cmd += ["-t", "4", *extra]
    if lowmem:
        cmd.append("--lowmem")
    proc = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=600,
                          stdin=subprocess.DEVNULL)
    assert proc.returncode == expect_rc, proc.stdout + proc.stderr
    return proc


def _collect(outdir):
    """Outputs worth pinning, normalized to be run-independent."""
    out = {}
    vcf = os.path.join(outdir, "debreak.vcf")
    if os.path.exists(vcf):
        with open(vcf) as f:
            out["debreak.vcf"] = _normalize_vcf(f.read())
    with open(os.path.join(outdir, "debreak-allsv-merged-final")) as f:
        out["debreak-allsv-merged-final"] = _sorted_lines(f.read())
    raw_dir = os.path.join(outdir, "sv_raw_calls")
    raw = "".join(open(os.path.join(raw_dir, n)).read() for n in sorted(os.listdir(raw_dir)))
    out["sv_raw_calls.txt"] = _sorted_lines(raw)
    with open(os.path.join(outdir, "log.txt")) as f:
        out["log.txt"] = f.read()
    return out


def _check_golden(name, outputs, update, required=True):
    gdir = os.path.join(GOLDEN, name)
    if update:
        shutil.rmtree(gdir, ignore_errors=True)
        os.makedirs(gdir)
        for fname, content in outputs.items():
            with open(os.path.join(gdir, fname), "w") as f:
                f.write(content)
        pytest.skip("golden files updated: " + gdir)
    if not os.path.isdir(gdir) and not required:
        pytest.skip(f"no golden files for '{name}' yet: generate them on Linux with "
                    "the external tools (CI: run workflow with update_golden=true)")
    assert os.path.isdir(gdir), "no golden files, run: pytest --update-golden"
    for fname, content in outputs.items():
        with open(os.path.join(gdir, fname)) as f:
            expected = f.read()
        assert content == expected, f"{name}/{fname} differs from golden"


SCENARIOS = {
    "default": [],
    "no_genotype": ["--no_genotype"],
    "tumor": ["--tumor"],
    "min_support_10": ["-m", "10"],
    "depth_30_maxcov_80": ["-d", "30", "--maxcov", "80"],
    "min_size_500": ["--min_size", "500"],
}


@BOTH_MODES
@pytest.mark.parametrize("name", SCENARIOS)
def test_pipeline_matches_golden(name, lowmem, tmp_path, update_golden):
    _run(tmp_path / "out", *SCENARIOS[name], lowmem=lowmem)
    _check_golden(name, _collect(tmp_path / "out"), update_golden and not lowmem)


@BOTH_MODES
def test_default_run_finds_all_six_simulated_svs(lowmem, tmp_path):
    """Sanity check independent of snapshots: the test BAM has 6 embedded SVs."""
    _run(tmp_path / "out", lowmem=lowmem)
    with open(tmp_path / "out" / "debreak.vcf") as f:
        records = [l.split("\t") for l in f if not l.startswith("#")]
    svtypes = sorted(r[7].split("SVTYPE=")[1].split(";")[0] for r in records)
    assert svtypes == ["DEL", "DEL", "DUP", "INS", "INS", "INV"]


def test_missing_input_exits_with_error(tmp_path):
    proc = subprocess.run(
        [sys.executable, LAUNCHER, "-o", str(tmp_path / "out")],
        cwd=REPO, capture_output=True, text=True,
    )
    assert proc.returncode == 1
    assert "No input file given" in proc.stdout


def test_rescue_dup_without_ref_logs_warning(tmp_path):
    _run(tmp_path / "out", "--rescue_dup")
    log = (tmp_path / "out" / "log.txt").read_text()
    assert "Reference genome is required for duplication rescue" in log


# --- modules needing external binaries (Linux / bioconda env) --------------

def _have(*tools):
    return all(shutil.which(t) for t in tools)


@BOTH_MODES
@pytest.mark.skipif(not _have("minimap2"), reason="needs minimap2")
def test_rescue_dup_matches_golden(lowmem, tmp_path, update_golden):
    _run(tmp_path / "out", "--rescue_dup", "-r", REF, lowmem=lowmem)
    _check_golden("rescue_dup", _collect(tmp_path / "out"), update_golden and not lowmem,
                  required=False)


@BOTH_MODES
@pytest.mark.skipif(not _have("minimap2", "wtdbg2", "wtpoa-cns", "bsalign"),
                    reason="needs minimap2, wtdbg2, bsalign")
def test_full_pipeline_matches_golden(lowmem, tmp_path, update_golden):
    _run(tmp_path / "out", "--poa", "--rescue_large_ins", "--rescue_dup", "-r", REF,
         lowmem=lowmem)
    _check_golden("full", _collect(tmp_path / "out"), update_golden and not lowmem,
                  required=False)


# --- translocations: synthetic 2-chromosome dataset (tests/make_tra_testdata.py)

TRA_DIR = os.path.join(REPO, "testdata", "tra")
TRA_SCENARIOS = {
    "tra_default": [],
    "tra_no_genotype": ["--no_genotype"],
    "tra_min_quality_61": ["--min_quality", "61"],  # all reads have MAPQ 60
    "tra_maxcov_10": ["--maxcov", "10"],
}


def _vcf_records(path):
    with open(path) as f:
        return [l.rstrip("\n").split("\t") for l in f if not l.startswith("#")]


def _info(rec):
    return dict(kv.split("=", 1) for kv in rec[7].split(";") if "=" in kv)


@BOTH_MODES
@pytest.mark.parametrize("name", TRA_SCENARIOS)
def test_tra_pipeline_matches_golden(name, lowmem, tmp_path, update_golden):
    _run(tmp_path / "out", *TRA_SCENARIOS[name], bam="tra.bam", cwd=TRA_DIR, lowmem=lowmem)
    _check_golden(name, _collect(tmp_path / "out"), update_golden and not lowmem)


@BOTH_MODES
def test_tra_both_translocations_detected(lowmem, tmp_path):
    _run(tmp_path / "out", bam="tra.bam", cwd=TRA_DIR, lowmem=lowmem)
    recs = _vcf_records(tmp_path / "out" / "debreak.vcf")
    assert [r[4] for r in recs] == ["<TRA>", "<TRA>"]
    # merged calls carry both breakpoints; positions from make_tra_testdata.EVENTS (+1)
    merged = (tmp_path / "out" / "translocation-merged").read_text().splitlines()
    assert sorted(tuple(l.split("\t")[:4]) for l in merged) == [
        ("chr1", "25001", "chr2", "30001"),
        ("chr1", "40001", "chr2", "10001"),
    ]
    for r in recs:
        info = _info(r)
        assert info["CHR2"] == "chr2"
        assert info["SUPPREAD"] == "25"
        assert r[9] == "1/1"


def test_tra_vcf_end_when_second_bp_before_pos(tmp_path):
    # SUSPECTED BUG (debreak_writevcf.py: `newrec.stop=int(sv[3])`): the TRA
    # breakpoint on CHR2 is stored as the record stop on the SAME contig.
    # When it is smaller than POS, htslib rejects it ("INFO/END=0 is smaller
    # than POS") and END silently becomes POS instead of the chr2 position.
    _run(tmp_path / "out", bam="tra.bam", cwd=TRA_DIR)
    recs = {r[1]: _info(r) for r in _vcf_records(tmp_path / "out" / "debreak.vcf")}
    assert recs["25002"]["END"] == "30001"  # bp2 > POS: correct
    assert recs["40002"]["END"] == "40002"  # bp2 = 10001 < POS: lost


def test_tra_testdata_is_reproducible(tmp_path):
    """The committed dataset matches what the generator produces."""
    import make_tra_testdata
    import pysam

    bam, fa = make_tra_testdata.generate(str(tmp_path))
    assert open(fa).read() == open(os.path.join(TRA_DIR, "tra_ref.fa")).read()
    with pysam.AlignmentFile(bam) as new, \
            pysam.AlignmentFile(os.path.join(TRA_DIR, "tra.bam")) as old:
        assert [a.to_string() for a in new] == [a.to_string() for a in old]


# --- --samlist input (SAM files instead of a sorted BAM) --------------------

@pytest.fixture(scope="module")
def samlist_dir(tmp_path_factory):
    """test_read.bam split into 2 name-grouped SAM files + samlist.txt."""
    import pysam

    d = tmp_path_factory.mktemp("samlist")
    with pysam.AlignmentFile(os.path.join(REPO, BAM), "rb") as bam:
        reads = sorted(bam.fetch(until_eof=True), key=lambda a: (a.query_name, a.flag))
        parts = [pysam.AlignmentFile(str(d / f"part{i}.sam"), "w", template=bam)
                 for i in range(2)]
        for a in reads:  # read names look like S1_123
            parts[int(a.query_name.split("_")[1]) % 2].write(a)
        for p in parts:
            p.close()
    (d / "samlist.txt").write_text("part0.sam\npart1.sam\n")
    return d


@BOTH_MODES
def test_samlist_matches_golden(samlist_dir, lowmem, tmp_path, update_golden):
    # BUG: --samlist always crashes when writing the VCF: writevcf_pysam()
    # gets the samlist *directory* as `bampath` and opens it with pysam.
    # Detection, clustering and debreak-allsv-merged-final still happen.
    proc = _run(tmp_path / "out", "--samlist", "samlist.txt", bam=None,
                cwd=samlist_dir, expect_rc=1, lowmem=lowmem)
    assert "Could not open alignment file" in proc.stderr
    assert not (tmp_path / "out" / "debreak.vcf").exists()
    _check_golden("samlist", _collect(tmp_path / "out"), update_golden and not lowmem)


def test_samlist_from_other_cwd_fails_on_depth_estimation(samlist_dir, tmp_path):
    # BUG: depth estimation opens samlist[0] relative to the CWD, not to the
    # samlist directory, so --samlist only works when run from that directory
    # (or when -d is given).
    proc = _run(tmp_path / "out", "--samlist", str(samlist_dir / "samlist.txt"),
                bam=None, expect_rc=1)
    assert "part0.sam" in proc.stderr
    log = (tmp_path / "out" / "log.txt").read_text()
    assert "Done detecting raw SV signals." in log
    assert "Estimated Sequencing Depth" not in log


def test_samlist_from_other_cwd_with_depth_reaches_vcf_step(samlist_dir, tmp_path):
    _run(tmp_path / "out", "--samlist", str(samlist_dir / "samlist.txt"), "-d", "50",
         bam=None, expect_rc=1)
    log = (tmp_path / "out" / "log.txt").read_text()
    assert "Warning:  No sorted bam file provided. Abort genotyping." in log
    assert log.rstrip().endswith("Start writing VCF file...")


def test_rescue_large_ins_requires_bam(samlist_dir, tmp_path):
    _run(tmp_path / "out", "--samlist", "samlist.txt", "--rescue_large_ins", "-r", "x.fa",
         bam=None, cwd=samlist_dir, expect_rc=1)
    log = (tmp_path / "out" / "log.txt").read_text()
    assert "Rescue_large_INS module requires sorted BAM input. Aborted." in log


# --- warning / option-handling branches of the `debreak` script ------------

def test_poa_without_ref_crashes_in_vcf_writer(tmp_path):
    # BUG: POA is skipped with a warning, but writevcf_pysam() is still called
    # with ifpoa=True and reffile=False; open(False) reads stdin (fd 0), the
    # reference dict stays empty and the first DEL raises KeyError.
    proc = _run(tmp_path / "out", "--poa", expect_rc=1)
    assert "KeyError: 'chr1'" in proc.stderr
    log = (tmp_path / "out" / "log.txt").read_text()
    assert "Warning: No reference genome provided. Abort breakpoint refinement." in log


@pytest.mark.parametrize("flag, warning", [
    ("--rescue_large_ins", "Warning: No reference provided. Rescue_large_INS module aborted."),
])
def test_modules_needing_ref_are_skipped_without_ref(flag, warning, tmp_path):
    _run(tmp_path / "out", flag)
    out = tmp_path / "out"
    assert warning in (out / "log.txt").read_text()
    # module skipped -> same SVs as the default run
    golden = os.path.join(GOLDEN, "default", "debreak.vcf")
    assert [r[:5] for r in _vcf_records(out / "debreak.vcf")] == \
        [r[:5] for r in _vcf_records(golden)]


def test_default_outpath_is_debreak_out_in_cwd(tmp_path):
    _run(None, bam=os.path.join(REPO, BAM), cwd=tmp_path)
    assert (tmp_path / "debreak-out" / "debreak.vcf").exists()


def test_outpath_without_trailing_slash(tmp_path):
    out = tmp_path / "nested" / "out"  # also: missing parents are created
    _run(None, "-o", str(out))
    assert (out / "debreak.vcf").exists()
    assert not (tmp_path / "nested" / "outdebreak.vcf").exists()


@BOTH_MODES
def test_skip_detect_reuses_existing_raw_calls(lowmem, tmp_path):
    out = tmp_path / "out"
    _run(out, lowmem=lowmem)
    first = (out / "debreak.vcf").read_text()
    (out / "map_depth" / "maplength_allchrom_debreak").unlink()
    _run(out, "--skip_detect", lowmem=lowmem)
    assert _normalize_vcf((out / "debreak.vcf").read_text()) == _normalize_vcf(first)
