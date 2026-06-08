#!/usr/bin/env python3

"""
find_softclip_sites_from_bam.py

samtools view SRR10403849_vs_UCSC_Hcal_v1.bam | ~/git/splice-leader/extra_scripts/find_softclip_sites_from_bam.py - > SRR10403849_vs_UCSC_Hcal_v1.softclips.tab
"""

import sys
import argparse
import re
from collections import Counter

cigar_re = re.compile(r"(\d+)([MIDNSHP=X])")


def is_reverse(flag):
    return flag & 16 != 0


def ref_consumed_by_cigar(cigar):
    total = 0
    for length, op in cigar_re.findall(cigar):
        if op in "MDN=X":
            total += int(length)
    return total


def fiveprime_softclip_site(rname, pos1, flag, cigar, min_clip):
    """
    Return BED-style site:
        chrom, start0, end0, strand

    SAM POS is 1-based leftmost aligned position.
    BED start is 0-based.
    """

    if cigar == "*":
        return None

    parts = cigar_re.findall(cigar)
    if not parts:
        return None

    strand = "-" if is_reverse(flag) else "+"

    left_len, left_op = int(parts[0][0]), parts[0][1]
    right_len, right_op = int(parts[-1][0]), parts[-1][1]

    ref_len = ref_consumed_by_cigar(cigar)

    # reference_start in 0-based BED coordinates
    ref_start0 = pos1 - 1

    # reference_end in 0-based half-open coordinates
    ref_end0 = ref_start0 + ref_len

    if not is_reverse(flag):
        if left_op == "S" and left_len >= min_clip:
            return rname, ref_start0, ref_start0 + 1, "+"
    else:
        if right_op == "S" and right_len >= min_clip:
            return rname, ref_end0 - 1, ref_end0, "-"

    return None


def main():
    ap = argparse.ArgumentParser(description="Extract 5' soft-clipped read starts from SAM stdin.")
    ap.add_argument('input_file', type = argparse.FileType('r'), default = '-', help="SAM file, - for stdin from BAM")
    ap.add_argument("-m", "--min-softclip", type=int, default=10, help="minimum length of the softclip in bp [10]")
    ap.add_argument("-c", "--min-count", type=int, default=1, help="minimum count of softclips at that base [1]")
    ap.add_argument("--bedgraph", action="store_true")
    args = ap.parse_args()

    counts = Counter()
    line_counter = 0
    total_softclips = 0

    sys.stderr.write("# Reading {}, tracking S of at least {}bp for count of at least {} \n".format(args.input_file.name, args.min_softclip, args.min_count ) )
    for line in args.input_file:
        line_counter += 1
        if not line or line.startswith("@"):
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) < 11:
            continue

        flag = int(fields[1])

        # skip unmapped reads
        if flag & 4:
            continue

        rname = fields[2]
        pos1 = int(fields[3])
        cigar = fields[5]

        site = fiveprime_softclip_site(
            rname=rname,
            pos1=pos1,
            flag=flag,
            cigar=cigar,
            min_clip=args.min_softclip,
        )

        if site is not None:
            total_softclips += 1
            counts[site] += 1

    for (chrom, start, end, strand), count in sorted(counts.items()):
        if count < args.min_count:
            continue

        if args.bedgraph:
            print(f"{chrom}\t{start}\t{end}\t{count}")
        else:
            name = f"softclip5p_count_{count}"
            print(f"{chrom}\t{start}\t{end}\t{name}\t{count}\t{strand}")

    sys.stderr.write("# Counted {} lines for {} leader sequences \n".format(line_counter, total_softclips ) )

if __name__ == "__main__":
    main()
