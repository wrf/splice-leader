#!/usr/bin/env python3

"""find_fused_cds_islands.py  v1  2026-06-08

~/git/splice-leader/extra_scripts/find_fused_cds_islands.py -g Hcv1av93_v10.leader_removed.gff -o Hcv1av93_v10.leader_removed.fused_cds.tab

Genes parsed: 14626
Transcripts parsed: 19660
Suspicious genes reported: 196
Output written to: Hcv1av93_v10.leader_removed.fused_cds.tab

"""

import argparse
import re
from collections import defaultdict


def parse_attrs(attr_text):
    attrs = {}
    for part in attr_text.strip().split(";"):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            attrs[k] = v
    return attrs


def parse_gff(path):
    genes = {}
    transcripts = {}
    cds_by_tx = defaultdict(list)
    tx_by_gene = defaultdict(list)

    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue

            seqid, source, ftype, start, end, score, strand, phase, attrs_raw = fields
            start, end = int(start), int(end)
            attrs = parse_attrs(attrs_raw)

            if ftype == "gene":
                gid = attrs.get("ID")
                if gid:
                    genes[gid] = {
                        "seqid": seqid,
                        "start": start,
                        "end": end,
                        "strand": strand,
                    }

            elif ftype in {"mRNA", "transcript"}:
                tid = attrs.get("ID")
                parents = attrs.get("Parent", "").split(",")
                if tid and parents:
                    transcripts[tid] = {
                        "seqid": seqid,
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "parents": parents,
                    }
                    for gid in parents:
                        tx_by_gene[gid].append(tid)

            elif ftype == "CDS":
                parents = attrs.get("Parent", "").split(",")
                for tid in parents:
                    if tid:
                        cds_by_tx[tid].append((seqid, start, end, strand))

    return genes, transcripts, tx_by_gene, cds_by_tx


def merge_intervals(intervals, max_gap=0):
    """
    Merge intervals if they overlap or are separated by <= max_gap.
    Intervals are (start, end).
    """
    if not intervals:
        return []

    intervals = sorted((min(s, e), max(s, e)) for s, e in intervals)
    merged = [intervals[0]]

    for s, e in intervals[1:]:
        last_s, last_e = merged[-1]
        if s <= last_e + max_gap + 1:
            merged[-1] = (last_s, max(last_e, e))
        else:
            merged.append((s, e))

    return merged


def intervals_overlap(a, b):
    """Return True if two closed intervals overlap."""
    a_start, a_end = a
    b_start, b_end = b
    return a_start <= b_end and b_start <= a_end


def find_suspicious_genes(
    genes,
    transcripts,
    tx_by_gene,
    cds_by_tx,
    min_islands=2,
):
    results = []

    for gid, txs in tx_by_gene.items():
        if gid not in genes:
            continue

        gene = genes[gid]

        tx_spans = {}

        for tid in txs:
            cds = cds_by_tx.get(tid, [])

            cds = [
                (seqid, start, end, strand)
                for seqid, start, end, strand in cds
                if seqid == gene["seqid"] and strand == gene["strand"]
            ]

            if not cds:
                continue

            cds_start = min(min(start, end) for _, start, end, _ in cds)
            cds_end = max(max(start, end) for _, start, end, _ in cds)

            tx_spans[tid] = (cds_start, cds_end)

        if len(tx_spans) < 2:
            continue

        # Build connected groups of transcript CDS spans.
        # Transcripts belong to the same group if their CDS spans overlap.
        groups = []

        for tid, span in sorted(tx_spans.items(), key=lambda x: x[1][0]):
            placed = False

            for group in groups:
                if any(intervals_overlap(span, other_span) for _, other_span in group):
                    group.append((tid, span))
                    placed = True
                    break

            if not placed:
                groups.append([(tid, span)])

        if len(groups) >= min_islands:
            results.append({
                "gene_id": gid,
                "seqid": gene["seqid"],
                "strand": gene["strand"],
                "gene_start": gene["start"],
                "gene_end": gene["end"],
                "num_transcripts": len(txs),
                "num_cds_span_groups": len(groups),
                "groups": groups,
                "transcript_spans": tx_spans,
            })

    return results


def write_report(results, outpath):
    with open(outpath, "w") as out:
        header = [
            "gene_id",
            "seqid",
            "strand",
            "gene_start",
            "gene_end",
            "num_transcripts",
            "num_cds_span_groups",
            "groups",
            "transcript_cds_spans",
        ]
        out.write("\t".join(header) + "\n")

        for r in results:
            group_strings = []

            for i, group in enumerate(r["groups"], start=1):
                members = []
                for tid, span in group:
                    members.append(f"{tid}:{span[0]}-{span[1]}")
                group_strings.append(f"group{i}=" + ",".join(members))

            tx_span_strings = [
                f"{tid}:{span[0]}-{span[1]}"
                for tid, span in sorted(r["transcript_spans"].items())
            ]

            row = [
                r["gene_id"],
                r["seqid"],
                r["strand"],
                str(r["gene_start"]),
                str(r["gene_end"]),
                str(r["num_transcripts"]),
                str(r["num_cds_span_groups"]),
                ";".join(group_strings),
                ";".join(tx_span_strings),
            ]

            out.write("\t".join(row) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Find genes with multiple non-overlapping CDS islands, often caused by fused genes with bridging transcripts."
    )
    parser.add_argument("-g", "--gff", required=True, help="Input GFF3 file")
    parser.add_argument("-o", "--output", required=True, help="Output TSV report")
    parser.add_argument(
        "--cds-merge-gap",
        type=int,
        default=1000,
        help="Merge CDS intervals into the same island if separated by this many bp or less. Default: 1000",
    )
    parser.add_argument(
        "--min-islands",
        type=int,
        default=2,
        help="Minimum number of CDS islands required to report a gene. Default: 2",
    )
    parser.add_argument(
        "--no-require-bridge",
        action="store_true",
        help="Report genes with multiple CDS islands even if no transcript spans first-to-last island.",
    )

    args = parser.parse_args()

    genes, transcripts, tx_by_gene, cds_by_tx = parse_gff(args.gff)

    results = find_suspicious_genes(
        genes,
        transcripts,
        tx_by_gene,
        cds_by_tx,
        min_islands=args.min_islands
    )

    write_report(results, args.output)

    print(f"Genes parsed: {len(genes)}")
    print(f"Transcripts parsed: {len(transcripts)}")
    print(f"Suspicious genes reported: {len(results)}")
    print(f"Output written to: {args.output}")


if __name__ == "__main__":
    main()
