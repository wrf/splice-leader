#!/usr/bin/env python3
# # v1.1 fix long first exon bug 2026-05-14

"""
detect_fused_gene_splits.py  last modified 2026-05-14

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

~/git/splice-leader/extra_scripts/detect_fused_gene_splits.py --gff Hcv1av93_v10.gff --genome Hcv1av93_v10_src.fa --motifs hcal_leader_motifs.fa --summary Hcv1av93_v10_splits_summary.tsv --splits-table Hcv1av93_v10_proposed_clusters.tsv

~/git/splice-leader/extra_scripts/detect_fused_gene_splits.py --gff GCF_026151205.1_MBARI_Bmic_1.0_genomic.gff --genome GCF_026151205.1_MBARI_Bmic_1.0_genomic.fna --motifs leader_seq.fasta --summary GCF_026151205.1_MBARI_Bmic_1.0_splits_summary.tsv --splits-table GCF_026151205.1_MBARI_Bmic_1.0_proposed_clusters.tsv

~/git/splice-leader/extra_scripts/detect_fused_gene_splits.py --gff GCA_048537945.1_crg_Mlei_v2_genomic.gff --genome GCA_048537945.1_crg_Mlei_v2_genomic.fna --motifs Mlei_leader_sequences.txt --summary GCA_048537945.1_crg_Mlei_v2_splits_summary.tsv --splits-table GCA_048537945.1_crg_Mlei_v2_proposed_clusters.tsv
"""

# ~/git/splice-leader/extra_scripts/detect_fused_gene_splits.py --gff Hcv1av93.gff --genome UCSC_Hcal_v1.fa --motifs hcal_leader_motifs.fa --summary Hcv1av93_splits_summary.tsv --splits-table Hcv1av93_proposed_clusters.tsv

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
	"""
	Match motif against the first exon sequence in transcript orientation.

	Removal rule should be based on whether the motif is at/near the
	transcript-oriented 3' end of the first exon.

	Returns:
	  motif_id, motif_len, match_len, identity, match_type,
	  match_start, match_end, upstream_bases, downstream_bases

	Coordinates are 0-based within the transcript-oriented exon sequence.
	"""

	seq = seq.upper()
	best_hit = None

	for motif_id, motif in motifs:
		motif = motif.upper()

		if not seq or not motif:
			continue

		motif_len = len(motif)
		exon_len = len(seq)

		# -------------------------------------------------------------
		# Case 1: exact motif occurs somewhere inside exon
		# -------------------------------------------------------------
		pos = seq.find(motif)
		if pos != -1:
			match_start = pos
			match_end = pos + motif_len
			downstream_bases = exon_len - match_end
			upstream_bases = match_start

			if downstream_bases == 0:
				match_type = "exon_3prime_end"
			else:
				match_type = "internal"

			hit = (
				motif_id,
				motif_len,
				motif_len,
				1.0,
				match_type,
				match_start,
				match_end,
				upstream_bases,
				downstream_bases,
			)

			if best_hit is None or downstream_bases < best_hit[8]:
				best_hit = hit

		# -------------------------------------------------------------
		# Case 2: exon is shorter than motif, but exon matches the
		# 3' end of the motif.
		# -------------------------------------------------------------
		L = min(exon_len, motif_len)

		exon_suffix = seq[-L:]
		motif_suffix = motif[-L:]

		matches = sum(
			1 for a, b in zip(exon_suffix, motif_suffix)
			if a == b
		)
		identity = matches / L

		if identity >= min_identity:
			match_start = exon_len - L
			match_end = exon_len
			downstream_bases = 0
			upstream_bases = match_start

			hit = (
				motif_id,
				motif_len,
				L,
				identity,
				"exon_3prime_end",
				match_start,
				match_end,
				upstream_bases,
				downstream_bases,
			)

			if best_hit is None or identity > best_hit[3]:
				best_hit = hit

	return best_hit

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
	max_extra_bases_for_removal,
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
				(
					motif_id,
					motif_len,
					match_len,
					identity,
					match_type,
					match_start,
					match_end,
					upstream_bases,
					downstream_bases,
				) = hit

				exon_len = fe["end"] - fe["start"] + 1

				remove_first_exon = downstream_bases < max_extra_bases_for_removal

				motif_hit_by_tx[tx_id] = {
					"motif_id": motif_id,
					"motif_len": motif_len,
					"match_len": match_len,
					"exon_len": exon_len,
					"identity": identity,
					"match_type": match_type,
					"match_start": match_start,
					"match_end": match_end,
					"upstream_bases": upstream_bases,
					"downstream_bases": downstream_bases,
					"remove_first_exon": remove_first_exon,
				}

				if remove_first_exon:
					false_first_exon_by_tx[tx_id] = fe

		# No motif-matching first exon in this gene
		if not motif_hit_by_tx:
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
		for tx_id, hit_info in motif_hit_by_tx.items():
			fe = first_exon(exons_by_tx[tx_id], transcripts[tx_id]["strand"])
			intervals = intervals_by_tx[tx_id]

			summary_rows.append(
				{
					"gene_id": gene_id,
					"transcript_id": tx_id,
					"motif_hit": hit_info["motif_id"],
					"false_first_exon": (
						f"{fe['seqid']}:{fe['start']}-{fe['end']}:{fe['strand']}"
					),
					"body_start": min(s for s, e in intervals),
					"body_end": max(e for s, e in intervals),
					"n_gene_transcripts": len(tx_ids),
					"n_motif_hit_transcripts": len(motif_hit_by_tx),
					"n_removable_first_exons": len(false_first_exon_by_tx),
					"n_clusters": len(clusters),
					"proposed_split": "yes" if is_proposed_split else "no",
					"cluster": cluster_lookup.get(tx_id, ""),
					"match_len": hit_info["match_len"],
					"match_start": hit_info["match_start"],
					"match_end": hit_info["match_end"],
					"upstream_bases": hit_info["upstream_bases"],
					"downstream_bases": hit_info["downstream_bases"],
					"exon_len": hit_info["exon_len"],
					"identity": f"{hit_info['identity']:.4f}",
					"match_type": hit_info["match_type"],
					"remove_first_exon": "yes" if hit_info["remove_first_exon"] else "no",
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
	columns = [
		"gene_id",
		"transcript_id",
		"motif_hit",
		"false_first_exon",
		"body_start",
		"body_end",
		"n_gene_transcripts",
		"n_motif_hit_transcripts",
		"n_removable_first_exons",
		"n_clusters",
		"proposed_split",
		"cluster",
		"motif_len",
		"match_len",
		"exon_len",
		"identity",
		"match_type",
		"match_start",
		"match_end",
		"upstream_bases",
		"downstream_bases",
		"remove_first_exon",
	]

	with open(path, "w") as out:
		out.write("\t".join(columns) + "\n")

		for row in summary_rows:
			out.write(
				"\t".join(str(row.get(col, "")) for col in columns)
				+ "\n"
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
	ap.add_argument("--max-extra-bases-for-removal", type=int, default=20, help="Only remove first exon if exon length is motif length plus at most this many extra bases. Default: 20" )
	ap.add_argument( "--min-identity", type=float, default=0.90, help="Minimum identity for motif match. Default: 0.90" )
	ap.add_argument( "--min-hit-txs", type=int, default=2, help="Minimum motif-hit transcripts required to flag a gene. Default: 2" )
	ap.add_argument( "--min-clusters", type=int, default=2, help="Minimum proposed transcript clusters required. Default: 2" )
	ap.add_argument( "--max-gap", type=int, default=10000, help="Maximum gap for body exons/transcripts to be grouped together. Default: 10000" )
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
		max_extra_bases_for_removal=args.max_extra_bases_for_removal,
	)

	write_summary(args.summary, summary_rows)
	write_splits_table(args.splits_table, results)

	print(f"Detected suspect fused genes: {len(results)}")
	print(f"Wrote summary: {args.summary}")
	print(f"Wrote splits table: {args.splits_table}")


if __name__ == "__main__":
	main()

example_results="""
Reading motifs from hcal_leader_motifs.fa
Found 380 motifs and reverse complements
Reading GFF from Hcv1av93.gff
Found 14591 genes and 20076 transcripts
Found motifs for 598 transcripts
Found 7 fused genes 
Detected suspect fused genes: 7
Wrote summary: Hcv1av93_splits_summary.tsv
Wrote splits table: Hcv1av93_proposed_clusters.tsv
Reading motifs from hcal_leader_motifs.fa
Found 380 motifs and reverse complements
Reading GFF from Hcv1av93_v10.gff
Found 14619 genes and 19660 transcripts
Found motifs for 601 transcripts
Found 7 fused genes 
Detected suspect fused genes: 7
Wrote summary: Hcv1av93_v10_splits_summary.tsv
Wrote splits table: Hcv1av93_v10_proposed_clusters.tsv


Reading motifs from leader_seq.fasta
Found 66 motifs and reverse complements
Reading GFF from GCF_026151205.1_MBARI_Bmic_1.0_genomic.gff
Found 15888 genes and 22067 transcripts
Found motifs for 1885 transcripts
Found 39 fused genes 
Detected suspect fused genes: 39
Wrote summary: GCF_026151205.1_MBARI_Bmic_1.0_splits_summary.tsv
Wrote splits table: GCF_026151205.1_MBARI_Bmic_1.0_proposed_clusters.tsv

Reading motifs from Mlei_leader_sequences.txt
Found 20 motifs and reverse complements
Reading GFF from GCA_048537945.1_crg_Mlei_v2_genomic.gff
Found 19625 genes and 26504 transcripts
Found motifs for 225 transcripts
Found 0 fused genes 
Detected suspect fused genes: 0
Wrote summary: GCA_048537945.1_crg_Mlei_v2_splits_summary.tsv
Wrote splits table: GCA_048537945.1_crg_Mlei_v2_proposed_clusters.tsv
"""
