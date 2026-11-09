"""Unit tests for clustering (debreak_merge_cluster), translocation merging
(debreak_merge_contig) and genotyping (debreak_genotype).

Characterization tests: they pin the current behaviour. Quirks are marked
`SUSPECTED QUIRK` / `SUSPECTED BUG`.
"""
import pysam
import pytest

import debreak_genotype as gt
import debreak_merge_cluster as mc
import debreak_merge_contig as mctg


def _ev(pos, size, i, mapq=60, kind="D-cigar", chrom="chr1"):
    """A raw SV signal line as produced by debreak_detect."""
    return f"{chrom}\t{pos}\t{size}\t{kind}\tr{i}_cigar1\t0\t{mapq}"


def _names(n, start=0):
    return ";".join(f"r{i}_cigar1" for i in range(start, start + n))


# ------------------------------------------------- mergeinfo_insertion_oneevent

def test_oneevent_summarizes_cluster():
    candi = [_ev(1000 + i, 100, i) for i in range(4)]
    assert mc.mergeinfo_insertion_oneevent(candi, 3) == [
        f"chr1\t1001\t100\t4\t4\t60.0\t{_names(4)}"
    ]


def test_oneevent_drops_too_long_outlier():
    candi = [_ev(1000, 100, i) for i in range(4)] + [_ev(1000, 300, 4)]
    (rec,) = mc.mergeinfo_insertion_oneevent(candi, 3)
    assert rec.split("\t")[2:5] == ["100", "4", "4"]
    assert "r4_cigar1" not in rec


def test_oneevent_drops_too_short_outlier():
    candi = [_ev(1000, 20, 0)] + [_ev(1000, 100, i) for i in range(1, 4)]
    (rec,) = mc.mergeinfo_insertion_oneevent(candi, 3)
    assert rec.split("\t")[2:4] == ["100", "3"]


def test_oneevent_below_min_support_gives_nothing():
    assert mc.mergeinfo_insertion_oneevent([_ev(1000, 100, i) for i in range(2)], 3) == []


def test_oneevent_needs_at_least_two_reads_even_if_min_support_is_one():
    assert mc.mergeinfo_insertion_oneevent([_ev(1000, 100, 0)], 1) == []


def test_oneevent_quality_is_floored_mean():
    candi = [_ev(1000, 100, 0, mapq=60), _ev(1000, 100, 1, mapq=59), _ev(1000, 100, 2, mapq=59)]
    assert mc.mergeinfo_insertion_oneevent(candi, 3)[0].split("\t")[5] == "59.0"


# ------------------------------------------------- mergeinfo_insertion / assign

def test_mergeinfo_single_allele_is_unique():
    candi = [_ev(1000, 100, i) for i in range(4)]
    (rec,) = mc.mergeinfo_insertion(candi, 3, "Deletion")
    assert rec.endswith("\tUnique\tDeletion")


def test_mergeinfo_two_alleles_are_compound():
    candi = [_ev(1000, 100, i) for i in range(4)] + [_ev(1000, 300, i) for i in range(4, 8)]
    recs = mc.mergeinfo_insertion(candi, 3, "Insertion")
    assert [r.split("\t")[2] for r in recs] == ["300", "100"]
    assert all(r.endswith("\tCompoundSV\tInsertion") for r in recs)


def test_mergeinfo_minor_allele_below_support_leaves_unique_major():
    candi = [_ev(1000, 100, i) for i in range(4)] + [_ev(1000, 300, i) for i in range(4, 6)]
    (rec,) = mc.mergeinfo_insertion(candi, 3, "Insertion")
    assert rec.split("\t")[2] == "100"
    assert rec.endswith("\tUnique\tInsertion")


def test_mergeinfo_below_support_gives_nothing():
    assert mc.mergeinfo_insertion([_ev(1000, 100, i) for i in range(2)], 3, "Deletion") == []


def test_assign_candi_insertion_splits_by_nearest_mean():
    candi = [_ev(1000, s, i) for i, s in enumerate([90, 110, 290, 310])]
    g1, g2, m1, m2 = mc.assign_candi_insertion(candi, 300, 100)
    assert [c.split("\t")[2] for c in g1] == ["290", "310"]
    assert [c.split("\t")[2] for c in g2] == ["90", "110"]
    assert (m1, m2) == (300, 100)


# ------------------------------------------------- cluster (DEL/DUP/INV)

def _cluster(tmp_path, events, svtype="del", mins=3, maxdepth=50, minmapq=None):
    mc.cluster(str(tmp_path) + "/", events, "chr1", 30_000, mins, maxdepth, svtype, minmapq)
    out = tmp_path / f"{svtype}-merged-cluster-chr1"
    return out.read_text().splitlines() if out.exists() else []


def test_cluster_empty_input_writes_nothing(tmp_path):
    assert mc.cluster(str(tmp_path) + "/", [], "chr1", 1000, 3, 50, "del", None) == 0
    assert list(tmp_path.iterdir()) == []


def test_cluster_small_deletions(tmp_path):
    recs = _cluster(tmp_path, [_ev(1000 + i, 200, i) for i in range(5)])
    assert len(recs) == 1
    f = recs[0].split("\t")
    assert (f[0], f[2], f[3]) == ("chr1", "200", "5")
    assert f[-2:] == ["Unique", "Deletion"]


@pytest.mark.parametrize("svtype, label", [("dup", "Duplication"), ("inv", "Inversion")])
def test_cluster_labels_by_svtype(tmp_path, svtype, label):
    recs = _cluster(tmp_path, [_ev(1000 + i, 200, i) for i in range(5)], svtype=svtype)
    assert recs[0].endswith("\t" + label)


def test_cluster_large_deletions_use_window_path(tmp_path):
    # > 3000bp: only the large-SV path (1600bp window) sees them
    recs = _cluster(tmp_path, [_ev(10_000 + i, 5000, i) for i in range(5)])
    assert recs == [f"chr1\t10002\t5000\t5\t5\t60.0\t{_names(5)}\tUnique\tDeletion"]


def test_cluster_large_and_small_path_deduplicated(tmp_path):
    # 2000 < size <= 3000: found by both paths, the large-path call is dropped
    recs = _cluster(tmp_path, [_ev(10_000 + i, 2500, i) for i in range(5)])
    assert len(recs) == 1


def test_cluster_two_separate_large_events(tmp_path):
    events = [_ev(5_000 + i, 5000, i) for i in range(4)] + \
             [_ev(20_000 + i, 4000, i) for i in range(4, 8)]
    recs = _cluster(tmp_path, events)
    assert [r.split("\t")[1:3] for r in recs] == [["5001", "5000"], ["20005", "4000"]]


@pytest.mark.parametrize("n, found", [(3, False), (4, True)])
def test_cluster_needs_more_than_three_overlapping_signals(tmp_path, n, found):
    # a region starts only where > 3 signals overlap (threshold=3), even if
    # min_support is lower
    recs = _cluster(tmp_path, [_ev(1000 + i, 200, i) for i in range(n)], mins=2)
    assert bool(recs) is found


@pytest.mark.parametrize("n, found", [(3, False), (4, True)])
def test_cluster_ins_needs_more_than_three_overlapping_signals(tmp_path, n, found):
    recs = _cluster_ins(tmp_path, [_ev(5000 + i, 150, i, kind="I-cigar") for i in range(n)], mins=2)
    assert bool(recs) is found


def test_cluster_region_deeper_than_maxdepth_is_dropped(tmp_path):
    assert _cluster(tmp_path, [_ev(1000 + i, 200, i) for i in range(5)], maxdepth=4) == []


def test_cluster_below_min_support_is_dropped(tmp_path):
    assert _cluster(tmp_path, [_ev(1000 + i, 200, i) for i in range(5)], mins=6) == []


@pytest.mark.parametrize("mapq, kept", [(60, True), (10, True), (5, False), (0, True)])
def test_cluster_default_mapq_filter(tmp_path, mapq, kept):
    # SUSPECTED QUIRK: without --min_quality, calls with mean MAPQ in [1, 10)
    # are dropped, but MAPQ 0 calls are kept.
    recs = _cluster(tmp_path, [_ev(1000 + i, 200, i, mapq=mapq) for i in range(5)])
    assert bool(recs) is kept


@pytest.mark.parametrize("mapq, kept", [(60, True), (20, True), (19, False), (0, False)])
def test_cluster_explicit_min_quality(tmp_path, mapq, kept):
    recs = _cluster(tmp_path, [_ev(1000 + i, 200, i, mapq=mapq) for i in range(5)], minmapq=20)
    assert bool(recs) is kept


# ------------------------------------------------- cluster_ins

def _cluster_ins(tmp_path, events, mins=3, maxdepth=50, minmapq=None):
    mc.cluster_ins(str(tmp_path) + "/", events, "chr1", 30_000, mins, maxdepth, minmapq)
    out = tmp_path / "ins-merged-cluster-chr1"
    return out.read_text().splitlines() if out.exists() else []


def test_cluster_ins_empty_input(tmp_path):
    assert mc.cluster_ins(str(tmp_path) + "/", [], "chr1", 1000, 3, 50, None) == 0


def test_cluster_ins_small_insertions_carry_region_columns(tmp_path):
    recs = _cluster_ins(tmp_path, [_ev(5000 + i, 150, i, kind="I-cigar") for i in range(5)])
    assert len(recs) == 1
    f = recs[0].split("\t")
    assert f[2:5] == ["150", "5", "5"]
    assert f[7:9] == ["Unique", "Insertion"]
    assert f[9] == "chr1"  # + region start, end, depth
    assert len(f) == 13


def test_cluster_ins_large_insertions_window_path(tmp_path):
    recs = _cluster_ins(tmp_path, [_ev(5000 + i, 5000, i, kind="I-cigar") for i in range(5)])
    assert recs == [f"chr1\t5002\t5000\t5\t5\t60.0\t{_names(5)}\tUnique\tInsertion"]


def test_cluster_ins_large_and_small_path_deduplicated(tmp_path):
    recs = _cluster_ins(tmp_path, [_ev(5000 + i, 2500, i, kind="I-cigar") for i in range(5)])
    assert len(recs) == 1
    assert len(recs[0].split("\t")) == 13  # the small-path record survives


def test_cluster_ins_appends_to_existing_file(tmp_path):
    # SUSPECTED QUIRK: unlike cluster(), cluster_ins() opens its output in
    # append mode, so re-running into the same directory duplicates calls.
    events = [_ev(5000 + i, 150, i, kind="I-cigar") for i in range(5)]
    _cluster_ins(tmp_path, list(events))
    assert len(_cluster_ins(tmp_path, list(events))) == 2


# ------------------------------------------------- translocation merging

def _tra(pos1, pos2, count, names, chr2="chr2", qual=60.0):
    return f"chr1\t{pos1}\t{chr2}\t{pos2}\t{count}\t{qual}\t0.0\t0.0\t{names}\tTranslocation"


def _raw_tra(pos1, pos2, i, chr2="chr2", mapq=60):
    """Raw TRA signal as built in `debreak`: chrom, pos, chr2, pos2, mapq, name."""
    return f"chr1\t{pos1}\t{chr2}\t{pos2}\t{mapq}\tr{i}_seg1"


def test_m_samechr_translocation_merges_nearby_calls():
    merged = mctg.m_samechr_translocation([_tra(1000, 5000, 3, "a"), _tra(1500, 5400, 1, "b")])
    assert merged == ["chr1\t1125\tchr2\t5100\t4\t60.0\t0.0\t0.0\ta;b\tTranslocation"]


def test_m_samechr_translocation_keeps_stronger_when_pos1_far():
    merged = mctg.m_samechr_translocation([_tra(1000, 5000, 3, "a"), _tra(9000, 5400, 5, "b")])
    assert merged == [_tra(9000, 5400, 5, "b")]


def test_m_samechr_translocation_keeps_first_when_stronger():
    merged = mctg.m_samechr_translocation([_tra(1000, 5000, 5, "a"), _tra(9000, 5400, 3, "b")])
    assert merged == [_tra(1000, 5000, 5, "a")]


def test_counttime_translocation_separate_windows():
    raw = [_raw_tra(1000 + i, 5000, i) for i in range(3)] + \
          [_raw_tra(9000 + i, 7000, i) for i in range(3, 6)] + \
          [_raw_tra(20000, 9000, 6)]  # lone call: below support
    recs = mctg.counttime_translocation(raw, 3)
    assert [r.split("\t")[1:4] for r in recs] == [["1001", "chr2", "5000"],
                                                  ["9004", "chr2", "7000"]]


def test_m_samechr_translocation_leaves_distinct_calls():
    calls = [_tra(1000, 5000, 3, "a"), _tra(1000, 9000, 3, "b"), _tra(1000, 5000, 3, "c", chr2="chr3")]
    assert sorted(mctg.m_samechr_translocation(list(calls))) == sorted(calls)


def test_mergeinfo_translocation_odd_count_median_and_outlier():
    candi = [_raw_tra(1000, 5000, 0), _raw_tra(1002, 5002, 1), _raw_tra(1004, 5004, 2),
             _raw_tra(1006, 9000, 3), _raw_tra(1008, 5006, 4)]
    (rec,) = mctg.mergeinfo_translocation(candi, 3)
    f = rec.split("\t")
    assert f[:5] == ["chr1", "1003", "chr2", "5003", "4"]  # pos2=9000 dropped
    assert "r3_seg1" not in f[8]


@pytest.mark.parametrize("offset, count", [(800, "5"), (801, "4")])
def test_mergeinfo_translocation_outlier_boundary_is_800(offset, count):
    candi = [_raw_tra(1000, 5000, i) for i in range(4)] + [_raw_tra(1000, 5000 + offset, 4)]
    (rec,) = mctg.mergeinfo_translocation(candi, 3)
    assert rec.split("\t")[4] == count


def test_mergeinfo_translocation_picks_most_common_chr2():
    candi = [_raw_tra(1000, 5000, i) for i in range(3)] + [_raw_tra(1000, 5000, 9, chr2="chr3")]
    (rec,) = mctg.mergeinfo_translocation(candi, 3)
    assert rec.split("\t")[2] == "chr2"


def test_mergeinfo_translocation_below_support():
    candi = [_raw_tra(1000, 5000, 0), _raw_tra(1000, 9000, 1), _raw_tra(1000, 20000, 2)]
    assert mctg.mergeinfo_translocation(candi, 3) == []


def test_merge_translocation_writes_file_and_applies_upper_bound(tmp_path):
    raw = [_raw_tra(1000 + i, 5000 + i, i) for i in range(4)]
    res = mctg.merge_translocation(3, str(tmp_path) + "/", list(raw), "chr1", True, None)
    assert len(res) == 1
    assert (tmp_path / "tra-merged-chr1").read_text() == res[0] + "\n"
    # upper bound: support must be <= 30 * min_support
    raw = [_raw_tra(1000 + i % 50, 5000, i) for i in range(100)]
    assert mctg.merge_translocation(3, str(tmp_path) + "/", raw, "chr1", True, None) == []


def test_merge_translocation_min_quality_filters_file_not_result(tmp_path):
    # SUSPECTED QUIRK: min_qual only filters what is written to the file;
    # the returned list still contains the low-quality call.
    raw = [_raw_tra(1000 + i, 5000, i, mapq=20) for i in range(4)]
    res = mctg.merge_translocation(3, str(tmp_path) + "/", raw, "chr1", True, 30)
    assert len(res) == 1
    assert (tmp_path / "tra-merged-chr1").read_text() == ""


def test_merge_translocation_empty_input(tmp_path):
    assert mctg.merge_translocation(3, str(tmp_path) + "/", [], "chr1", True, None) is True


# ------------------------------------------------- genotyping
# Tiny BAMs with exactly controlled coverage. chr1 is 20 kbp; SV at 5000.

def _make_bam(path, reads, chroms=(("chr1", 20_000), ("chr2", 20_000))):
    """reads: list of (chrom, start, cigar, flag)."""
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": c, "LN": l} for c, l in chroms]})
    unsorted = str(path) + ".unsorted.bam"
    with pysam.AlignmentFile(unsorted, "wb", header=header) as f:
        for i, (chrom, start, cigar, flag) in enumerate(reads):
            a = pysam.AlignedSegment(header)
            a.query_name, a.flag, a.reference_name = f"q{i}", flag, chrom
            a.reference_start, a.mapping_quality, a.cigarstring = start, 60, cigar
            a.query_sequence = "A" * a.infer_query_length()
            f.write(a)
    pysam.sort("-o", str(path), unsorted)
    pysam.index(str(path))
    return str(path)


def _spanning(n, start=4000, chrom="chr1", cigar="2000M"):
    return [(chrom, start, cigar, 0)] * n


def _genotype(tmp_path, func, sv_lines, reads, highcov=50, chrom="chr1"):
    bam = _make_bam(tmp_path / "x.bam", reads)
    calls = tmp_path / "calls"
    calls.write_text("".join(l + "\n" for l in sv_lines))
    ret = func(bam, str(calls), chrom, highcov)
    out = tmp_path / ("calls-gt" if chrom == "all" else f"calls-gt-{chrom}")
    return ret, (out.read_text().splitlines() if out.exists() else None)


def _del(supp, pos=5000, size=200, chrom="chr1"):
    return f"{chrom}\t{pos}\t{size}\t{supp}\t{supp}\t60.0\tnames\tUnique\tDeletion"


@pytest.mark.parametrize("supp, expected", [(6, "GT=1/1"), (5, "GT=1/0"), (1, "GT=1/0")])
def test_genotype_del_threshold_is_60_percent_of_local_cov(tmp_path, supp, expected):
    _, out = _genotype(tmp_path, gt.genotype_del, [_del(supp)], _spanning(10))
    assert out == [_del(supp) + "\t" + expected]


def test_genotype_del_support_above_highcov_dropped(tmp_path):
    _, out = _genotype(tmp_path, gt.genotype_del, [_del(51)], _spanning(10))
    assert out == []


@pytest.mark.parametrize("start", [4000, 5250])  # spans left side+right / right side only
def test_genotype_del_flank_coverage_above_2x_highcov_dropped(tmp_path, start):
    _, out = _genotype(tmp_path, gt.genotype_del, [_del(3)], _spanning(10, start=start), highcov=4)
    assert out == []


def test_genotype_del_local_cov_between_highcov_and_2x_dropped(tmp_path):
    _, out = _genotype(tmp_path, gt.genotype_del, [_del(3)], _spanning(10), highcov=6)
    assert out == []


def test_genotype_del_no_calls_for_chrom(tmp_path):
    ret, out = _genotype(tmp_path, gt.genotype_del, [_del(3, chrom="chr2")], _spanning(10))
    assert ret == 0 and out is None


def test_genotype_del_all_chroms(tmp_path):
    _, out = _genotype(tmp_path, gt.genotype_del, [_del(6)], _spanning(10), chrom="all")
    assert out == [_del(6) + "\tGT=1/1"]


def test_genotype_filter_del(tmp_path):
    _, out = _genotype(tmp_path, gt.genotype_filter_del, [_del(1)], _spanning(10))
    assert out == [_del(1) + "\tPASS"]


@pytest.mark.parametrize("start", [4000, 5250])
def test_genotype_filter_del_high_flank_cov_dropped(tmp_path, start):
    _, out = _genotype(tmp_path, gt.genotype_filter_del, [_del(1)], _spanning(10, start=start),
                       highcov=4)
    assert out == []


def test_genotype_filter_del_no_calls_for_chrom(tmp_path):
    ret, out = _genotype(tmp_path, gt.genotype_filter_del, [_del(1, chrom="chr2")], _spanning(10))
    assert ret == 0 and out is None


def _ins(supp, pos=5000, size=300):
    return f"chr1\t{pos}\t{size}\t{supp}\t{supp}\t60.0\tnames\tUnique\tInsertion"


def test_genotype_ins_direct_1_1(tmp_path):
    _, out = _genotype(tmp_path, gt.genotype_ins, [_ins(6)], _spanning(10))
    assert out == [_ins(6) + "\tGT=1/1"]


def test_genotype_ins_rescan_without_insertions_gives_1_0(tmp_path):
    _, out = _genotype(tmp_path, gt.genotype_ins, [_ins(3)], _spanning(10))
    assert out == [_ins(3) + "\tGT=1/0"]


def test_genotype_ins_rescan_counts_cigar_insertions(tmp_path):
    # numsupp=3 < 0.6*16, but rescanning finds 12/16 reads with a ~300bp INS
    reads = _spanning(12, cigar="1000M300I1000M") + _spanning(4)
    _, out = _genotype(tmp_path, gt.genotype_ins, [_ins(3)], reads)
    assert out == [_ins(3) + "\tGT=1/1"]


def test_genotype_ins_rescan_ignores_reverse_strand_reads(tmp_path):
    # SUSPECTED QUIRK: the rescan skips every read with flag != 0, so
    # reverse-strand (flag 16) reads never count as support.
    reads = [("chr1", 4000, "1000M300I1000M", 16)] * 12 + _spanning(4)
    _, out = _genotype(tmp_path, gt.genotype_ins, [_ins(3)], reads)
    assert out == [_ins(3) + "\tGT=1/0"]


def test_genotype_ins_large_insertion_counts_left_clipped_reads(tmp_path):
    # svsize >= 500: reads starting at the breakpoint with a soft clip count (6/8)
    reads = [("chr1", 5000, "500S2000M", 0)] * 6 + _spanning(2)
    _, out = _genotype(tmp_path, gt.genotype_ins, [_ins(2, size=800)], reads)
    assert out == [_ins(2, size=800) + "\tGT=1/1"]


def test_genotype_ins_large_insertion_ignores_right_clipped_reads(tmp_path):
    # SUSPECTED BUG (operator precedence in genotype_ins):
    #   `if svsize<500 and abs(ref_start-start)<=100 or abs(start-ref_end)<=100: continue`
    # parses as `(A and B) or C`, so reads ENDING near the breakpoint are always
    # skipped, even for large insertions. The right-clip support branch below it
    # is therefore unreachable here: 6/8 supporting reads still give 1/0.
    reads = [("chr1", 3000, "2000M500S", 0)] * 6 + _spanning(2)
    _, out = _genotype(tmp_path, gt.genotype_ins, [_ins(2, size=800)], reads)
    assert out == [_ins(2, size=800) + "\tGT=1/0"]


@pytest.mark.parametrize("highcov, supp", [(50, 51), (4, 3)])
def test_genotype_ins_dropped(tmp_path, highcov, supp):
    _, out = _genotype(tmp_path, gt.genotype_ins, [_ins(supp)], _spanning(10), highcov=highcov)
    assert out == []


def test_genotype_ins_right_flank_too_deep_dropped(tmp_path):
    _, out = _genotype(tmp_path, gt.genotype_ins, [_ins(3)], _spanning(10, start=5100), highcov=4)
    assert out == []


def test_genotype_filter_ins(tmp_path):
    _, out = _genotype(tmp_path, gt.genotype_filter_ins, [_ins(1)], _spanning(10))
    assert out == [_ins(1) + "\tPASS"]


@pytest.mark.parametrize("start", [4000, 5100])
def test_genotype_filter_ins_high_cov_dropped(tmp_path, start):
    _, out = _genotype(tmp_path, gt.genotype_filter_ins, [_ins(1)], _spanning(10, start=start),
                       highcov=4)
    assert out == []


def _tra_call(supp=6):
    return f"chr1\t5000\tchr2\t8000\t{supp}\t60.0\t0.0\t0.0\tnames\tTranslocation"


def _tra_reads(n1=10, n2=10):
    return _spanning(n1) + _spanning(n2, start=7000, chrom="chr2")


@pytest.mark.parametrize("supp, expected", [(6, "GT=1/1"), (5, "GT=1/0")])
def test_genotype_tra(tmp_path, supp, expected):
    _, out = _genotype(tmp_path, gt.genotype_tra, [_tra_call(supp)], _tra_reads())
    assert out == [_tra_call(supp) + "\t" + expected]


def test_genotype_tra_local_cov_is_max_over_both_breakpoints(tmp_path):
    _, out = _genotype(tmp_path, gt.genotype_tra, [_tra_call(6)], _tra_reads(n2=20))
    assert out == [_tra_call(6) + "\tGT=1/0"]  # 6 < 0.6 * 20


@pytest.mark.parametrize("reads", [
    _spanning(10),                                   # chr1 both flanks deep
    _spanning(10, start=5100),                       # chr1 right flank deep
    _spanning(10, start=7000, chrom="chr2"),         # chr2 both flanks deep
    _spanning(10, start=8100, chrom="chr2"),         # chr2 right flank deep
])
def test_genotype_tra_deep_flank_dropped(tmp_path, reads):
    _, out = _genotype(tmp_path, gt.genotype_tra, [_tra_call(2)], reads, highcov=4)
    assert out == []


def test_genotype_tra_support_above_highcov_dropped(tmp_path):
    _, out = _genotype(tmp_path, gt.genotype_tra, [_tra_call(51)], _tra_reads())
    assert out == []


def test_genotype_filter_tra(tmp_path):
    _, out = _genotype(tmp_path, gt.genotype_filter_tra, [_tra_call(1)], _tra_reads())
    assert out == [_tra_call(1) + "\tPASS"]


@pytest.mark.parametrize("reads", [
    _spanning(10),
    _spanning(10, start=5100),
    _spanning(10, start=7000, chrom="chr2"),
    _spanning(10, start=8100, chrom="chr2"),
])
def test_genotype_filter_tra_deep_flank_dropped(tmp_path, reads):
    _, out = _genotype(tmp_path, gt.genotype_filter_tra, [_tra_call(1)], reads, highcov=4)
    assert out == []
