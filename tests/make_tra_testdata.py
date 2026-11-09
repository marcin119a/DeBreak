"""Generate a small synthetic dataset with two translocations.

testdata/ only has one chromosome, so no translocation code path is ever
exercised. This script builds a deterministic 2-chromosome dataset:

  - reference: chr1, chr2 (50 kbp each, random sequence, fixed seed)
  - ~10x background coverage of perfectly matching 2 kbp reads
  - TRA_A: chr1:25000 joined to chr2:30000 (primary on chr1, SA on chr2)
  - TRA_B: chr2:10000 joined to chr1:40000 (primary on chr2, SA on chr1)

Each junction read is written as a primary alignment (soft-clipped) plus a
supplementary alignment (hard-clipped), both carrying SA tags, as minimap2 does.

Output is committed in testdata/tra/. Regenerate with:
    python tests/make_tra_testdata.py
"""
import os
import random
import sys

import pysam

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "testdata", "tra")

SEED = 20261109
CHROM_LEN = 50_000
READ_LEN = 2_000
BACKGROUND_DEPTH = 10
JUNCTION_READS = 25

# (name, left chrom, left breakpoint, right chrom, right breakpoint)
EVENTS = [
    ("traA", "chr1", 25_000, "chr2", 30_000),
    ("traB", "chr2", 10_000, "chr1", 40_000),
]


def _randseq(rng, n):
    return "".join("ACGT"[int(rng.random() * 4)] for _ in range(n))


def _segment(header, name, flag, chrom, start, cigar, seq, sa):
    a = pysam.AlignedSegment(header)
    a.query_name = name
    a.flag = flag
    a.reference_name = chrom
    a.reference_start = start
    a.mapping_quality = 60
    a.cigarstring = cigar
    a.query_sequence = seq
    a.set_tag("SA", sa)
    return a


def generate(out=OUT):
    rng = random.Random(SEED)
    ref = {c: _randseq(rng, CHROM_LEN) for c in ("chr1", "chr2")}
    os.makedirs(out, exist_ok=True)

    fa = os.path.join(out, "tra_ref.fa")
    with open(fa, "w") as f:
        for chrom, seq in ref.items():
            f.write(f">{chrom}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")
    pysam.faidx(fa)

    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": c, "LN": CHROM_LEN} for c in ref],
    })
    reads = []

    # background coverage
    n_bg = BACKGROUND_DEPTH * CHROM_LEN // READ_LEN
    for chrom in ref:
        for i in range(n_bg):
            start = int(rng.random() * (CHROM_LEN - READ_LEN))
            a = pysam.AlignedSegment(header)
            a.query_name = f"bg_{chrom}_{i}"
            a.flag = 0 if rng.random() < 0.5 else 16
            a.reference_name = chrom
            a.reference_start = start
            a.mapping_quality = 60
            a.cigarstring = f"{READ_LEN}M"
            a.query_sequence = ref[chrom][start:start + READ_LEN]
            reads.append(a)

    # junction reads: left part on chrom1 ending at bp1, right part on chrom2 from bp2
    for ev, chr1, bp1, chr2, bp2 in EVENTS:
        for i in range(JUNCTION_READS):
            left = 700 + int(rng.random() * 600)
            right = READ_LEN - left
            name = f"{ev}_{i}"
            seq = ref[chr1][bp1 - left:bp1] + ref[chr2][bp2:bp2 + right]
            p_start, s_start = bp1 - left, bp2
            # SA format: rname,pos(1-based),strand,CIGAR,mapQ,NM;
            sa_to_supp = f"{chr2},{s_start + 1},+,{left}S{right}M,60,0;"
            sa_to_prim = f"{chr1},{p_start + 1},+,{left}M{right}S,60,0;"
            reads.append(_segment(header, name, 0, chr1, p_start,
                                  f"{left}M{right}S", seq, sa_to_supp))
            reads.append(_segment(header, name, 2048, chr2, s_start,
                                  f"{left}H{right}M", seq[left:], sa_to_prim))

    unsorted = os.path.join(out, "tra_unsorted.bam")
    with pysam.AlignmentFile(unsorted, "wb", header=header) as f:
        for a in reads:
            f.write(a)
    bam = os.path.join(out, "tra.bam")
    pysam.sort("-o", bam, unsorted)
    os.remove(unsorted)
    pysam.index(bam)
    return bam, fa


if __name__ == "__main__":
    print(generate(sys.argv[1] if len(sys.argv) > 1 else OUT))
