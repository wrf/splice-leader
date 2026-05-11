#!/usr/bin/env python3

"""
detect_fused_gene_splits.py  last modified 2026-05-11

Detect genes with transcripts sharing a false first exon motif, then propose
transcript clusters for splitting the gene.

Inputs:
  --gff	   GFF3 annotation
  --genome	genome FASTA
  --motifs	FASTA of known false/shared first-exon motifs

Outputs:
  --summary   detailed TSV
  --splits	Python dictionary text suitable for hardcoded SPLITS

Requires:
  pip install biopython

~/git/splice-leader/extra_scripts/detect_fused_gene_splits.py --gff Hcv1av93.gff --genome UCSC_Hcal_v1.fa --motifs hcal_leader_motifs.fa --summary Hcv1av93_splits_summary.tsv --splits-table Hcv1av93_proposed_clusters.tsv
"""

import sys
import argparse
from collections import defaultdict
from Bio import SeqIO

##############################

def parse_attrs(attr_string):
	attrs = {}
	for part in attr_string.strip().split(";"):
		if not part:
			continue
		if "=" in part:
			k, v = part.split("=", 1)
			attrs[k] = v
	return attrs

##############################

def parent_ids(attrs):
	parent = attrs.get("Parent", "")
	if not parent:
		return []
	return parent.split(",")

##############################

def revcomp(seq):
	table = str.maketrans("ACGTacgt", "TGCAtgca")
	return seq.translate(table)[::-1]

##############################

def load_genome(path):
	return {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(path, "fasta")}

##############################

def load_motifs(path):
	motifs = []
	print(f"Reading motifs from {path}", file=sys.stderr )
	for rec in SeqIO.parse(path, "fasta"):
		seq = str(rec.seq).upper()
		motifs.append((rec.id, seq))
		motifs.append((rec.id + "_revcomp", revcomp(seq).upper()))
	print("Found {} motifs and reverse complements".format( len(motifs) ), file=sys.stderr )
	return motifs

##############################

def motif_matches(seq, motifs, min_identity=0.90):
	seq = seq.upper()

	for motif_id, motif in motifs:
		motif = motif.upper()

		if motif in seq or seq in motif:
			return motif_id

		L = min(len(seq), len(motif))
		if L == 0:
			continue

		matches = sum(1 for a, b in zip(seq[:L], motif[:L]) if a == b)
		identity = matches / L

		if identity >= min_identity:
			return motif_id

	return None

##############################

def read_gff(path):
	genes = {}
	transcripts = {}
	tx_by_gene = defaultdict(list)
	exons_by_tx = defaultdict(list)
	print(f"Reading GFF from {path}", file=sys.stderr )
	with open(path) as fh:
		for line in fh:
			if line.startswith("#") or not line.strip():
				continue

			parts = line.rstrip("\n").split("\t")
			if len(parts) != 9:
				continue

			seqid, source, ftype, start, end, score, strand, phase, attrs_s = parts
			start = int(start)
			end = int(end)
			attrs = parse_attrs(attrs_s)

			feature = {
				"seqid": seqid,
				"source": source,
				"type": ftype,
				"start": start,
				"end": end,
				"strand": strand,
				"attrs": attrs,
				"line": line.rstrip("\n"),
			}

			if ftype == "gene" and "ID" in attrs:
				genes[attrs["ID"]] = feature

			elif ftype in {"transcript", "mRNA"} and "ID" in attrs:
				tx_id = attrs["ID"]
				transcripts[tx_id] = feature

				for parent in parent_ids(attrs):
					tx_by_gene[parent].append(tx_id)

			elif ftype == "exon":
				for parent in parent_ids(attrs):
					exons_by_tx[parent].append(feature)
	print("Found {} genes and {} transcripts".format( len(genes), len(transcripts) ), file=sys.stderr )
	return genes, transcripts, tx_by_gene, exons_by_tx

##############################

def first_exon(exons, strand):
	if not exons:
		return None

	if strand == "-":
		return max(exons, key=lambda x: x["end"])
	else:
		return min(exons, key=lambda x: x["start"])

##############################

def get_exon_seq(genome, exon):
	seq = genome[exon["seqid"]][exon["start"] - 1:exon["end"]]
	if exon["strand"] == "-":
		seq = revcomp(seq)
	return seq.upper()

##############################

def interval_overlap_or_close(a, b, max_gap):
	a_start, a_end = a
	b_start, b_end = b

	if a_start <= b_end and b_start <= a_end:
		return True

	if b_start > a_end and b_start - a_end <= max_gap:
		return True

	if a_start > b_end and a_start - b_end <= max_gap:
		return True

	return False

##############################

def transcript_body_intervals(tx_id, transcripts, exons_by_tx, false_first_exon_by_tx):
	"""
	Use exon intervals, excluding the false first exon if present.
	If no exons remain, fall back to the transcript span.
	"""
	tx = transcripts[tx_id]
	intervals = []

	false_exon = false_first_exon_by_tx.get(tx_id)

	for exon in exons_by_tx.get(tx_id, []):
		if false_exon:
			same_exon = (
				exon["seqid"] == false_exon["seqid"]
				and exon["start"] == false_exon["start"]
				and exon["end"] == false_exon["end"]
				and exon["strand"] == false_exon["strand"]
			)
			if same_exon:
				continue

		intervals.append((exon["start"], exon["end"]))

	if not intervals:
		intervals = [(tx["start"], tx["end"])]

	return sorted(intervals)

##############################

def transcripts_connected(tx1, tx2, intervals_by_tx, max_gap):
	for a in intervals_by_tx[tx1]:
		for b in intervals_by_tx[tx2]:
			if interval_overlap_or_close(a, b, max_gap):
				return True
	return False

##############################

def cluster_transcripts(tx_ids, intervals_by_tx, max_gap):
	tx_ids = list(tx_ids)
	parent = {x: x for x in tx_ids}

	def find(x):
		while parent[x] != x:
			parent[x] = parent[parent[x]]
			x = parent[x]
		return x

	def union(a, b):
		ra = find(a)
		rb = find(b)
		if ra != rb:
			parent[rb] = ra

	for i in range(len(tx_ids)):
		for j in range(i + 1, len(tx_ids)):
			if transcripts_connected(tx_ids[i], tx_ids[j], intervals_by_tx, max_gap):
				union(tx_ids[i], tx_ids[j])

	clusters = defaultdict(list)
	for tx in tx_ids:
		clusters[find(tx)].append(tx)

	return list(clusters.values())

##############################

def cluster_span(cluster, intervals_by_tx):
	starts = []
	ends = []

	for tx_id in cluster:
		for s, e in intervals_by_tx[tx_id]:
			starts.append(s)
			ends.append(e)

	return min(starts), max(ends)

##############################

def detect_splits(
	genes,
	transcripts,
	tx_by_gene,
	exons_by_tx,
	genome,
	motifs,
	min_identity,
	min_hit_txs,
	min_clusters,
	max_gap,
):
	results = []
	summary_rows = []

	for gene_id, tx_ids in tx_by_gene.items():
		if gene_id not in genes:
			continue

		false_first_exon_by_tx = {}
		motif_hit_by_tx = {}

		# -------------------------------------------------------------
		# First pass: detect motif-matching first exons
		# -------------------------------------------------------------
		for tx_id in tx_ids:
			if tx_id not in transcripts:
				continue

			tx = transcripts[tx_id]
			exons = exons_by_tx.get(tx_id, [])

			fe = first_exon(exons, tx["strand"])
			if fe is None:
				continue

			seq = get_exon_seq(genome, fe)
			hit = motif_matches(seq, motifs, min_identity=min_identity)

			if hit:
				false_first_exon_by_tx[tx_id] = fe
				motif_hit_by_tx[tx_id] = hit

		# No motif-matching first exon in this gene
		if not false_first_exon_by_tx:
			continue

		# -------------------------------------------------------------
		# Cluster all transcripts in the gene after excluding the
		# suspicious first exon from motif-hit transcripts.
		# -------------------------------------------------------------
		intervals_by_tx = {
			tx_id: transcript_body_intervals(
				tx_id,
				transcripts,
				exons_by_tx,
				false_first_exon_by_tx,
			)
			for tx_id in tx_ids
			if tx_id in transcripts
		}

		clusters = cluster_transcripts(
			list(intervals_by_tx.keys()),
			intervals_by_tx,
			max_gap=max_gap,
		)

		clusters = sorted(
			clusters,
			key=lambda c: cluster_span(c, intervals_by_tx)[0],
		)

		is_proposed_split = (
			len(false_first_exon_by_tx) >= min_hit_txs
			and len(clusters) >= min_clusters
		)

		cluster_lookup = {}
		for cluster_i, cluster in enumerate(clusters, start=1):
			for tx_id in cluster:
				cluster_lookup[tx_id] = cluster_i

		# -------------------------------------------------------------
		# Add summary rows for every transcript with a motif-hit first exon,
		# even when the gene is not proposed for splitting.
		# -------------------------------------------------------------
		for tx_id, fe in false_first_exon_by_tx.items():
			intervals = intervals_by_tx[tx_id]

			summary_rows.append(
				{
					"gene_id": gene_id,
					"transcript_id": tx_id,
					"motif_hit": motif_hit_by_tx.get(tx_id, ""),
					"false_first_exon": (
						f"{fe['seqid']}:{fe['start']}-{fe['end']}:{fe['strand']}"
					),
					"body_start": min(s for s, e in intervals),
					"body_end": max(e for s, e in intervals),
					"n_gene_transcripts": len(tx_ids),
					"n_motif_hit_transcripts": len(false_first_exon_by_tx),
					"n_clusters": len(clusters),
					"proposed_split": "yes" if is_proposed_split else "no",
					"cluster": cluster_lookup.get(tx_id, ""),
				}
			)

		# -------------------------------------------------------------
		# Only genes passing fusion/split criteria go into the splits table.
		# -------------------------------------------------------------
		if is_proposed_split:
			results.append(
				{
					"gene_id": gene_id,
					"tx_ids": tx_ids,
					"false_first_exon_by_tx": false_first_exon_by_tx,
					"motif_hit_by_tx": motif_hit_by_tx,
					"intervals_by_tx": intervals_by_tx,
					"clusters": clusters,
				}
			)
	print("Found motifs for {} transcripts".format( len(summary_rows) ), file=sys.stderr )
	print("Found {} fused genes ".format( len(results) ), file=sys.stderr )
	return results, summary_rows

##############################

def write_summary(path, summary_rows):
	with open(path, "w") as out:
		out.write(
			"gene_id\ttranscript_id\tmotif_hit\tfalse_first_exon\t"
			"body_start\tbody_end\tn_gene_transcripts\t"
			"n_motif_hit_transcripts\tn_clusters\tproposed_split\tcluster\n"
		)

		for row in summary_rows:
			out.write(
				f"{row['gene_id']}\t"
				f"{row['transcript_id']}\t"
				f"{row['motif_hit']}\t"
				f"{row['false_first_exon']}\t"
				f"{row['body_start']}\t"
				f"{row['body_end']}\t"
				f"{row['n_gene_transcripts']}\t"
				f"{row['n_motif_hit_transcripts']}\t"
				f"{row['n_clusters']}\t"
				f"{row['proposed_split']}\t"
				f"{row['cluster']}\n"
			)

##############################

def write_splits_table(path, results):
	with open(path, "w") as out:
		out.write("old_gene_id\tcluster\ttranscript_id\n")

		for result in results:
			gene_id = result["gene_id"]
			intervals_by_tx = result["intervals_by_tx"]

			for cluster_i, cluster in enumerate(result["clusters"], start=1):

				cluster = sorted(
					cluster,
					key=lambda tx: intervals_by_tx[tx][0][0],
				)

				for tx_id in cluster:
					out.write(
						f"{gene_id}\t{cluster_i}\t{tx_id}\n"
					)

##############################

def main():
	ap = argparse.ArgumentParser()

	ap.add_argument("--gff", required=True)
	ap.add_argument("--genome", required=True)
	ap.add_argument("--motifs", required=True)

	ap.add_argument("--summary", required=True)
	ap.add_argument("--splits-table", required=True, help="Output TSV containing proposed transcript clusters")

	ap.add_argument(
		"--min-identity",
		type=float,
		default=0.90,
		help="Minimum identity for motif match. Default: 0.90",
	)

	ap.add_argument(
		"--min-hit-txs",
		type=int,
		default=2,
		help="Minimum motif-hit transcripts required to flag a gene. Default: 2",
	)

	ap.add_argument(
		"--min-clusters",
		type=int,
		default=2,
		help="Minimum proposed transcript clusters required. Default: 2",
	)

	ap.add_argument(
		"--max-gap",
		type=int,
		default=10000,
		help="Maximum gap for body exons/transcripts to be grouped together. Default: 10000",
	)

	args = ap.parse_args()

	genome = load_genome(args.genome)
	motifs = load_motifs(args.motifs)
	genes, transcripts, tx_by_gene, exons_by_tx = read_gff(args.gff)

	results, summary_rows = detect_splits(
		genes=genes,
		transcripts=transcripts,
		tx_by_gene=tx_by_gene,
		exons_by_tx=exons_by_tx,
		genome=genome,
		motifs=motifs,
		min_identity=args.min_identity,
		min_hit_txs=args.min_hit_txs,
		min_clusters=args.min_clusters,
		max_gap=args.max_gap,
	)

	write_summary(args.summary, summary_rows)
	write_splits_table(args.splits_table, results)

	print(f"Detected suspect fused genes: {len(results)}")
	print(f"Wrote summary: {args.summary}")
	print(f"Wrote splits table: {args.splits_table}")


if __name__ == "__main__":
	main()
