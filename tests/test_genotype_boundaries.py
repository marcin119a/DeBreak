"""Boundary tests for debreak_genotype, written to kill surviving mutants.

Every coverage window, threshold and column index is probed exactly at its
edge (value and value +/- 1), so that an off-by-one change in the code flips
the result. Reads are short (10 bp) so that each one touches a single window.

pysam's count(contig, start, stop) counts reads overlapping [start, stop).
"""
import pysam
import pytest

import debreak_genotype as gt

CHROMS = (("chr1", 20_000), ("chr2", 20_000))


def _bam(path, reads):
    """reads: iterable of (chrom, start, cigar[, flag])."""
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": c, "LN": l} for c, l in CHROMS]})
    unsorted = str(path) + ".u.bam"
    with pysam.AlignmentFile(unsorted, "wb", header=header) as f:
        for i, r in enumerate(reads):
            chrom, start, cigar, flag = (*r, 0) if len(r) == 3 else r
            a = pysam.AlignedSegment(header)
            a.query_name, a.flag, a.reference_name = f"q{i}", flag, chrom
            a.reference_start, a.mapping_quality, a.cigarstring = start, 60, cigar
            a.query_sequence = "A" * a.infer_query_length()
            f.write(a)
    pysam.sort("-o", str(path), unsorted)
    pysam.index(str(path))
    return str(path)


def _run(tmp_path, func, lines, reads, highcov, chrom="chr1"):
    bam = _bam(tmp_path / "x.bam", reads)
    calls = tmp_path / "calls"
    calls.write_text("".join(l + "\n" for l in lines))
    ret = func(bam, str(calls), chrom, highcov)
    out = tmp_path / ("calls-gt" if chrom == "all" else f"calls-gt-{chrom}")
    return ret, (out.read_text().splitlines() if out.exists() else None)


def _reads_ending_at(chrom, end, n=3):
    return [(chrom, end - 10, "10M")] * n


def _reads_starting_at(chrom, start, n=3):
    return [(chrom, start, "10M")] * n


def _edge_reads(chrom, lo, hi, edge, inside):
    """3 short reads just inside / just outside one edge of window [lo, hi)."""
    if edge == "lo":
        return _reads_ending_at(chrom, lo + 1 if inside else lo)
    return _reads_starting_at(chrom, hi - 1 if inside else hi)


# Records: column 3 (numsupp) deliberately differs from column 4.
def DEL(pos=5000, size=200, supp=1, chrom="chr1"):
    return f"{chrom}\t{pos}\t{size}\t{supp}\t99\t60.0\tn\tUnique\tDeletion"


def INS(pos=5000, size=300, supp=1, chrom="chr1"):
    return f"{chrom}\t{pos}\t{size}\t{supp}\t99\t60.0\tn\tUnique\tInsertion"


def TRA(bp1=5000, bp2=8000, supp=1, chrom="chr1", chr2="chr2"):
    return f"{chrom}\t{bp1}\t{chr2}\t{bp2}\t{supp}\t60.0\t0.0\t0.0\tn\tTranslocation"


# window edges per function: (chrom, lo, hi) for the default record
DEL_WINDOWS = [("chr1", 4850, 4950), ("chr1", 5250, 5350)]  # stop = 5200
INS_WINDOWS = [("chr1", 4850, 4950), ("chr1", 5050, 5150)]
TRA_WINDOWS = [("chr1", 4850, 4950), ("chr1", 5050, 5150),
               ("chr2", 7850, 7950), ("chr2", 8050, 8150)]

PASS_FUNCS = {  # function -> (record, windows, label written when kept)
    "genotype_filter_del": (DEL(), DEL_WINDOWS, "PASS"),
    "genotype_del": (DEL(), DEL_WINDOWS, "GT=1/1"),
    "genotype_filter_ins": (INS(), INS_WINDOWS, "PASS"),
    "genotype_ins": (INS(), INS_WINDOWS, "GT=1/1"),
    "genotype_filter_tra": (TRA(), TRA_WINDOWS, "PASS"),
    "genotype_tra": (TRA(), TRA_WINDOWS, "GT=1/1"),
}

EDGE_CASES = [
    (name, w, edge, inside)
    for name, (_, windows, _) in PASS_FUNCS.items()
    for w in range(len(windows))
    for edge in ("lo", "hi")
    for inside in (True, False)
]


@pytest.mark.parametrize("name, w, edge, inside", EDGE_CASES)
def test_coverage_window_edges(tmp_path, name, w, edge, inside):
    """3 reads inside a window -> local cov 3 > highcov 2 -> call dropped;
    the same reads 1bp outside the window -> cov 0 -> call kept."""
    rec, windows, label = PASS_FUNCS[name]
    reads = _edge_reads(*windows[w], edge, inside)
    _, out = _run(tmp_path, getattr(gt, name), [rec], reads, highcov=2)
    assert out == ([] if inside else [rec + "\t" + label])


@pytest.mark.parametrize("name", PASS_FUNCS)
def test_returns_true_and_keeps_every_record_in_order(tmp_path, name):
    rec, _, label = PASS_FUNCS[name]
    second = rec.replace("\t5000\t", "\t12000\t", 1)
    ret, out = _run(tmp_path, getattr(gt, name), [rec, second], [], highcov=2)
    assert ret is True
    assert out == [rec + "\t" + label, second + "\t" + label]


@pytest.mark.parametrize("name", PASS_FUNCS)
def test_filters_records_by_chromosome(tmp_path, name):
    rec, _, label = PASS_FUNCS[name]
    other = "chr2" + rec[4:]
    (tmp_path / "a").mkdir(), (tmp_path / "b").mkdir()
    ret, out = _run(tmp_path / "a", getattr(gt, name), [other, rec], [], highcov=2)
    assert out == [rec + "\t" + label]
    ret, out = _run(tmp_path / "b", getattr(gt, name), [other], [], highcov=2)
    assert ret == 0 and out is None  # nothing for chr1: no file at all


@pytest.mark.parametrize("name", PASS_FUNCS)
def test_all_mode_keeps_every_chromosome(tmp_path, name):
    rec, _, label = PASS_FUNCS[name]
    other = "chr2" + rec[4:] if "tra" not in name else rec  # TRA: chr2 is CHR2
    lines = [rec, other] if other != rec else [rec]
    ret, out = _run(tmp_path, getattr(gt, name), lines, [], highcov=2, chrom="all")
    assert ret is True
    assert out == [l + "\t" + label for l in lines]


@pytest.mark.parametrize("name", PASS_FUNCS)
def test_counts_only_reads_on_the_records_chromosome(tmp_path, name):
    # reads at the same coordinates on another chromosome must not count
    rec, windows, label = PASS_FUNCS[name]
    chrom, lo, _ = windows[0]
    other = "chr2" if chrom == "chr1" else "chr1"
    reads = _reads_starting_at(other, lo, n=5)
    if "tra" in name:  # chr2 is used by TRA itself; use a third-chrom-free check
        reads = _reads_starting_at("chr2", 4850, n=5)  # chr2:4850 is not a TRA window
    _, out = _run(tmp_path, getattr(gt, name), [rec], reads, highcov=2)
    assert out == [rec + "\t" + label]


@pytest.mark.parametrize("name", PASS_FUNCS)
@pytest.mark.parametrize("cov, kept", [(2, True), (3, False)])
def test_local_coverage_equal_to_highcov_is_kept(tmp_path, name, cov, kept):
    rec, windows, label = PASS_FUNCS[name]
    chrom, lo, _ = windows[0]
    reads = _reads_starting_at(chrom, lo, n=cov)
    _, out = _run(tmp_path, getattr(gt, name), [rec.replace("\t1\t", "\t2\t", 1)], reads,
                  highcov=2)
    assert bool(out) is kept


# ------------------------------------------------- support column / thresholds

@pytest.mark.parametrize("func, mk", [(gt.genotype_del, DEL), (gt.genotype_ins, INS)])
@pytest.mark.parametrize("supp, kept", [(2, True), (3, False)])
def test_support_equal_to_highcov_is_kept(tmp_path, func, mk, supp, kept):
    _, out = _run(tmp_path, func, [mk(supp=supp)], [], highcov=2)
    assert bool(out) is kept


@pytest.mark.parametrize("supp, kept", [(2, True), (3, False)])
def test_tra_support_equal_to_highcov_is_kept(tmp_path, supp, kept):
    _, out = _run(tmp_path, gt.genotype_tra, [TRA(supp=supp)], [], highcov=2)
    assert bool(out) is kept


def _cov(chrom, lo, n):
    return _reads_starting_at(chrom, lo, n=n)


@pytest.mark.parametrize("supp, gt_label", [(3, "GT=1/1"), (2, "GT=1/0")])
def test_del_genotype_uses_support_column_3(tmp_path, supp, gt_label):
    # local cov 5: 3 >= 0.6*5 -> 1/1, 2 < 3 -> 1/0 (column 4 is always 99)
    _, out = _run(tmp_path, gt.genotype_del, [DEL(supp=supp)], _cov("chr1", 4850, 5), highcov=10)
    assert out == [DEL(supp=supp) + "\t" + gt_label]


@pytest.mark.parametrize("supp, gt_label", [(3, "GT=1/1"), (2, "GT=1/0")])
def test_tra_genotype_uses_support_column_5(tmp_path, supp, gt_label):
    _, out = _run(tmp_path, gt.genotype_tra, [TRA(supp=supp)], _cov("chr2", 8050, 5), highcov=10)
    assert out == [TRA(supp=supp) + "\t" + gt_label]


def test_del_near_contig_start_clamps_window_to_zero(tmp_path):
    # SV at 100: left window is [max(0, -50), 50) -> a read on [0, 1) counts
    rec = DEL(pos=100)
    _, out = _run(tmp_path, gt.genotype_del, [rec], [("chr1", 0, "1M")] * 3, highcov=2)
    assert out == []


# ------------------------------------------------- genotype_ins: rescan phase
# INS at 5000, numsupp 1 < 0.6 * local cov, so the reads overlapping 5000 are
# rescanned. Spanning reads cover [4000, 6000) -> they also count in both
# coverage windows. highcov 50.

SPAN = ("chr1", 4000, "2000M")


def _ins_gt(tmp_path, reads, size=300, supp=1):
    _, out = _run(tmp_path, gt.genotype_ins, [INS(size=size, supp=supp)], reads, highcov=50)
    assert out and len(out) == 1
    return out[0].rsplit("\t", 1)[1]


def _ins_read(ins_len, ins_at=5000, start=4000, total=2000, pre_cigar=""):
    """Read [start, start+total) with an insertion of ins_len at ref pos ins_at."""
    left = ins_at - start
    return ("chr1", start, f"{pre_cigar}{left}M{ins_len}I{total - left}M")


@pytest.mark.parametrize("ins_len, counted", [
    (210, True), (209, False),   # 0.7 * 300 = 210
    (429, True), (430, False),   # 1.43 * 300 = 429
])
def test_ins_rescan_insertion_size_ratio_bounds(tmp_path, ins_len, counted):
    reads = [_ins_read(ins_len)] * 6 + [SPAN] * 4  # 6/10 support if counted
    assert _ins_gt(tmp_path, reads) == ("GT=1/1" if counted else "GT=1/0")


@pytest.mark.parametrize("ins_at, counted", [(5500, True), (5501, False), (4500, True)])
def test_ins_rescan_insertion_distance_bound_500(tmp_path, ins_at, counted):
    reads = [_ins_read(300, ins_at=ins_at, start=3800, total=2400)] * 6 + [SPAN] * 4
    assert _ins_gt(tmp_path, reads) == ("GT=1/1" if counted else "GT=1/0")


def test_ins_rescan_far_insertion_of_right_size_not_counted(tmp_path):
    reads = [_ins_read(300, ins_at=5800, start=3800, total=2400)] * 6 + [SPAN] * 4
    assert _ins_gt(tmp_path, reads) == "GT=1/0"


def test_ins_rescan_offset_counts_deletions_before_insertion(tmp_path):
    # 100M100D...: the D advances the reference offset, so the insertion lands
    # at 3800 + 100 + 100 + 1000 = 5000 (counted). Without the D it would be at 4900.
    reads = [("chr1", 3800, "100M100D1000M300I1000M")] * 6 + [SPAN] * 4
    assert _ins_gt(tmp_path, reads) == "GT=1/1"


def test_ins_rescan_offset_ignores_soft_clips(tmp_path):
    # leading soft clip does not move the reference offset
    reads = [("chr1", 4000, "300S1000M300I1000M")] * 6 + [SPAN] * 4
    assert _ins_gt(tmp_path, reads) == "GT=1/1"


def test_ins_rescan_offset_ignores_ref_skips(tmp_path):
    # SUSPECTED QUIRK: N (op 3) consumes reference but is not added to the
    # offset (only M=0 and D=2 are): the insertion really at 3400+600+1000=5000
    # is computed at 4400, more than 500bp away -> not counted.
    reads = [("chr1", 3400, "600N1000M300I1000M")] * 6 + [SPAN] * 4
    assert _ins_gt(tmp_path, reads) == "GT=1/0"


def test_ins_rescan_counts_each_read_once(tmp_path):
    # two matching insertions in one read still count as one supporting read:
    # 3/10 -> 1/0 (would be 6/10 -> 1/1 if counted twice)
    two = ("chr1", 4000, "1000M300I100M300I900M")
    reads = [two] * 3 + [SPAN] * 7
    assert _ins_gt(tmp_path, reads) == "GT=1/0"


@pytest.mark.parametrize("n_ins, label", [(6, "GT=1/1"), (5, "GT=1/0")])
def test_ins_rescan_ratio_threshold_60_percent(tmp_path, n_ins, label):
    reads = [_ins_read(300)] * n_ins + [SPAN] * (10 - n_ins)
    assert _ins_gt(tmp_path, reads) == label


def test_ins_rescan_counts_all_overlapping_reads(tmp_path):
    reads = [_ins_read(300)] * 3 + [SPAN] * 7  # 3/10
    assert _ins_gt(tmp_path, reads) == "GT=1/0"
    reads = [_ins_read(300)] * 3 + [SPAN] * 2  # 3/5
    assert _ins_gt(tmp_path, reads) == "GT=1/1"


@pytest.mark.parametrize("start, overlaps", [(5000, True), (5001, False)])
def test_ins_rescan_fetch_right_edge(tmp_path, start, overlaps):
    # large INS (no start-distance skip): a plain read starting at 5000
    # overlaps the fetch window [4999, 5001) and dilutes support 5/5 -> 5/10
    clipped = [("chr1", 5000, "500S1000M")] * 5
    plain = [("chr1", start, "1000M")] * 5
    label = _ins_gt(tmp_path, clipped + plain, size=800)
    assert label == ("GT=1/1" if not overlaps else "GT=1/1")  # 5/10 >= 0.5 for large INS
    # make dilution decisive: 4 clipped + 5 plain -> 4/9 < 0.5 only if plain overlap
    label = _ins_gt(tmp_path, clipped[:4] + plain, size=800)
    assert label == ("GT=1/0" if overlaps else "GT=1/1")


@pytest.mark.parametrize("end, overlaps", [(4999 + 1, True), (4999, False)])
def test_ins_rescan_fetch_left_edge(tmp_path, end, overlaps):
    # reads ending near 5000 are skipped anyway (precedence quirk), so the left
    # fetch edge is unobservable: both cases look the same.
    reads = [_ins_read(300)] * 6 + [("chr1", end - 1000, "1000M")] * 10
    assert _ins_gt(tmp_path, reads) == "GT=1/1"


@pytest.mark.parametrize("size, skipped", [(499, True), (500, False)])
def test_ins_rescan_start_distance_skip_only_for_small_ins(tmp_path, size, skipped):
    # plain reads starting 100bp after the INS are skipped only when size < 500
    near = [("chr1", 4900, "1000M")] * 6
    reads = [_ins_read(size)] * 4 + near
    label = _ins_gt(tmp_path, reads, size=size)
    # skipped: 4/4 -> 1/1; not skipped: 4/10 -> 1/0 (also < 0.5 for large)
    assert label == ("GT=1/1" if skipped else "GT=1/0")


@pytest.mark.parametrize("delta, skipped", [(100, True), (101, False), (-100, True), (-101, False)])
def test_ins_rescan_start_distance_bound_100(tmp_path, delta, skipped):
    near = [("chr1", 5000 + delta - 0, "1000M")] * 6 if delta < 0 else \
           [("chr1", 5000 - 1000 + 1, "1000M")] * 0
    # reads starting at 5000+delta (delta<0 -> before the INS, still overlapping it)
    near = [("chr1", 5000 + delta, "2000M")] * 6
    reads = [_ins_read(300)] * 4 + near
    label = _ins_gt(tmp_path, reads)
    if delta > 0:  # starting after 5000 -> not in fetch window at all
        assert label == "GT=1/1"
    else:
        assert label == ("GT=1/1" if skipped else "GT=1/0")


@pytest.mark.parametrize("delta, skipped", [(100, True), (101, False)])
def test_ins_rescan_end_distance_bound_100(tmp_path, delta, skipped):
    # reads ending 100/101bp after the INS start: skipped within 100bp
    near = [("chr1", 5000 + delta - 2000, "2000M")] * 6
    reads = [_ins_read(300)] * 4 + near
    assert _ins_gt(tmp_path, reads) == ("GT=1/1" if skipped else "GT=1/0")


def test_ins_rescan_ignores_reads_far_from_start(tmp_path):
    # SUSPECTED QUIRK: precedence `(A and B) or C` - for size>=500 a read ending
    # within 100bp is skipped (C), a read starting near the INS is not (A false)
    near_start = [("chr1", 5000, "1000M")] * 4  # A false for size 800 -> counted
    clipped = [("chr1", 5000, "500S1000M")] * 4
    assert _ins_gt(tmp_path, clipped + near_start, size=800) == "GT=1/1"  # 4/8 >= 0.5


@pytest.mark.parametrize("cigar, counted", [
    ("500S1000M", True), ("500H1000M", True), ("1000M", False), ("10I1000M", False),
])
def test_ins_rescan_left_clip_types(tmp_path, cigar, counted):
    reads = [("chr1", 5000, cigar)] * 6 + [SPAN] * 4
    assert _ins_gt(tmp_path, reads, size=800) == ("GT=1/1" if counted else "GT=1/0")


@pytest.mark.parametrize("start, counted", [(5250, True), (5251, False), (4750, True)])
def test_ins_rescan_left_clip_distance_bound_250(tmp_path, start, counted):
    # left-clipped reads must overlap 5000 to be fetched: use long reads that
    # start before 5000 when start < 5000; reads starting after 5000 are not fetched
    if start > 5000:
        pytest.skip("reads starting after the INS are outside the fetch window")
    reads = [("chr1", start, "500S1000M")] * 6 + [SPAN] * 4
    assert _ins_gt(tmp_path, reads, size=800) == ("GT=1/1" if counted else "GT=1/0")


@pytest.mark.parametrize("start, counted", [(4749, False), (4750, True)])
def test_ins_rescan_left_clip_distance_bound_250_before(tmp_path, start, counted):
    reads = [("chr1", start, "500S1000M")] * 6 + [SPAN] * 4
    assert _ins_gt(tmp_path, reads, size=800) == ("GT=1/1" if counted else "GT=1/0")


@pytest.mark.parametrize("end, counted", [(5101, True), (5250, True), (5251, False)])
def test_ins_rescan_right_clip_counts_beyond_100bp(tmp_path, end, counted):
    # reads ending 101..250bp after the INS with a right clip count as support
    # (those ending within 100bp are skipped, see the precedence quirk)
    reads = [("chr1", end - 1000, "1000M500S")] * 6 + [SPAN] * 4
    assert _ins_gt(tmp_path, reads, size=800) == ("GT=1/1" if counted else "GT=1/0")


@pytest.mark.parametrize("cigar, counted", [
    ("1000M500S", True), ("1000M500H", True), ("1000M", False), ("1000M10I", False),
])
def test_ins_rescan_right_clip_types(tmp_path, cigar, counted):
    reads = [("chr1", 4200, cigar)] * 6 + [SPAN] * 4  # ends at 5200
    assert _ins_gt(tmp_path, reads, size=800) == ("GT=1/1" if counted else "GT=1/0")


def test_ins_rescan_right_clip_needs_last_op(tmp_path):
    # only the LAST cigar op is checked: clip followed by M is not a right clip
    reads = [("chr1", 4200, "990M500S10M")] * 6 + [SPAN] * 4
    assert _ins_gt(tmp_path, reads, size=800) == "GT=1/0"


@pytest.mark.parametrize("size, n_clip, label", [
    (500, 5, "GT=1/1"),   # large: 5/10 >= 0.5
    (500, 4, "GT=1/0"),
    (800, 5, "GT=1/1"),
])
def test_ins_rescan_large_ins_threshold_50_percent(tmp_path, size, n_clip, label):
    reads = [("chr1", 5000, "500S1000M")] * n_clip + [SPAN] * (10 - n_clip)
    assert _ins_gt(tmp_path, reads, size=size) == label


def test_ins_rescan_small_ins_threshold_is_60_percent(tmp_path):
    # size 499: 5/10 is not enough (needs 0.6), clips are not counted anyway
    reads = [_ins_read(499)] * 5 + [SPAN] * 5
    assert _ins_gt(tmp_path, reads, size=499) == "GT=1/0"


def test_ins_rescan_only_forward_unflagged_reads(tmp_path):
    reads = [_ins_read(300) + (16,)] * 6 + [_ins_read(300) + (0,)] * 4 + [SPAN] * 6
    assert _ins_gt(tmp_path, reads) == "GT=1/0"  # only 4/10 flag-0 reads counted


def test_ins_local_cov_between_highcov_and_2x_writes_nothing(tmp_path):
    # localcov 3 with highcov 2: no rescan, no output (and no NameError)
    _, out = _run(tmp_path, gt.genotype_ins, [INS()], _cov("chr1", 4850, 3), highcov=2)
    assert out == []


# ------------------------------------------------- second round: remaining survivors

def _shift(rec, name, by=7000):
    """Same record moved `by` bp to the right (both breakpoints for TRA)."""
    f = rec.split("\t")
    f[1] = str(int(f[1]) + by)
    if "tra" in name:
        f[3] = str(int(f[3]) + by)
    return "\t".join(f)


WINDOW_CASES = [(name, w) for name, (_, ws, _) in PASS_FUNCS.items() for w in range(len(ws))]


@pytest.mark.parametrize("name, w", WINDOW_CASES)
def test_deep_window_skips_only_that_record(tmp_path, name, w):
    # 5 reads > 2*highcov in one window: the record is skipped (`continue`),
    # processing must go on with the next record
    rec, windows, label = PASS_FUNCS[name]
    chrom, lo, _ = windows[w]
    nxt = _shift(rec, name)
    _, out = _run(tmp_path, getattr(gt, name), [rec, nxt], _cov(chrom, lo, 5), highcov=2)
    assert out == [nxt + "\t" + label]


@pytest.mark.parametrize("func, mk", [(gt.genotype_del, DEL), (gt.genotype_ins, INS)])
def test_high_support_skips_only_that_record(tmp_path, func, mk):
    _, out = _run(tmp_path, func, [mk(supp=3), mk(pos=12000, supp=1)], [], highcov=2)
    assert out == [mk(pos=12000, supp=1) + "\tGT=1/1"]


def test_tra_high_support_skips_only_that_record(tmp_path):
    _, out = _run(tmp_path, gt.genotype_tra, [TRA(supp=3), TRA(bp1=12000, bp2=15000)], [],
                  highcov=2)
    assert out == [TRA(bp1=12000, bp2=15000) + "\tGT=1/1"]


@pytest.mark.parametrize("name, w", WINDOW_CASES)
@pytest.mark.parametrize("cov, kept", [(2, True), (3, False)])
def test_each_window_coverage_equal_to_highcov_is_kept(tmp_path, name, w, cov, kept):
    rec, windows, _ = PASS_FUNCS[name]
    chrom, lo, _ = windows[w]
    rec = rec.replace("\t1\t", "\t2\t", 1)  # support 2: genotype 1/1 at cov 2
    _, out = _run(tmp_path, getattr(gt, name), [rec], _cov(chrom, lo, cov), highcov=2)
    assert bool(out) is kept


def test_ins_rescan_offset_sums_all_ops_before_insertion(tmp_path):
    # offset = 500 + 10 + 490 = 1000 -> insertion at 4000 + 1000 = 5000
    reads = [("chr1", 4000, "500M10D490M300I1000M")] * 6 + [SPAN] * 4
    assert _ins_gt(tmp_path, reads) == "GT=1/1"


def test_ins_rescan_ignores_reads_left_of_the_insertion(tmp_path):
    # reads entirely before 4999 are outside the fetch window [4999, 5001)
    far_left = [("chr1", 1000, "1000M")] * 20
    reads = [_ins_read(300)] * 6 + [SPAN] * 4 + far_left
    assert _ins_gt(tmp_path, reads) == "GT=1/1"


def test_ins_rescan_runs_when_local_cov_equals_highcov(tmp_path):
    # local cov == local cov2 == highcov == 10, support 1 -> rescan -> 6/10 -> 1/1
    reads = [_ins_read(300)] * 6 + [SPAN] * 4
    _, out = _run(tmp_path, gt.genotype_ins, [INS(supp=1)], reads, highcov=10)
    assert out == [INS(supp=1) + "\tGT=1/1"]


def test_ins_rescan_scans_reads_after_a_clipped_one(tmp_path):
    # supporting clipped reads come first in coordinate order is not guaranteed;
    # mix them so that a `break` after the first counted read changes the result
    reads = [("chr1", 4990, "500S1000M")] * 3 + [("chr1", 5000, "500S1000M")] * 3 + [SPAN] * 4
    assert _ins_gt(tmp_path, reads, size=800) == "GT=1/1"  # 6/10


def test_ins_rescan_continues_after_a_flagged_read(tmp_path):
    # a reverse-strand read sorting BEFORE the others is skipped and the scan
    # goes on: 3/10 -> 1/0. (A `break` there would leave 0/0, and since
    # 0 >= 0.6*0 is true, it would report 1/1.)
    reads = [("chr1", 3900, "2100M", 16)] + [_ins_read(300)] * 3 + [SPAN] * 7
    assert _ins_gt(tmp_path, reads) == "GT=1/0"


def test_ins_rescan_left_fetch_edge_is_unobservable():
    """Documentation only: mutants moving the left fetch edge (start-1 -> start-2)
    are EQUIVALENT: the extra read would end at start-1, and reads ending within
    100bp of the INS are always skipped by the precedence quirk."""
