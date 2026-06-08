#!/usr/bin/env python3
# v1.1 fix long first exon bug 2026-05-14

"""
split_fused_gff_by_false_initial_exon.py  last modified 2026-05-14

~/git/splice-leader/extra_scripts/split_fused_gff_by_false_initial_exon.py --gff Hcv1av93_v10.gff --splits-table Hcv1av93_v10_proposed_clusters.tsv --motif-summary Hcv1av93_v10_splits_summary.tsv --out-gff Hcv1av93_v10.leader_removed.gff --repair-summary Hcv1av93_v10.leader_removed.repair_summary.tsv
"""

# ~/git/splice-leader/extra_scripts/split_fused_gff_by_false_initial_exon.py --gff Hcv1av93.gff --splits-table Hcv1av93_proposed_clusters.tsv --motif-summary Hcv1av93_splits_summary.tsv --out-gff Hcv1av93.leader_removed.gff --repair-summary Hcv1av93.leader_removed.repair_summary.tsv


import sys
import argparse
import string
from collections import defaultdict

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

def format_attrs(attrs):
	return ";".join(f"{k}={v}" for k, v in attrs.items())

##############################

def parent_ids(attrs):
	if "Parent" not in attrs:
		return []
	return attrs["Parent"].split(",")

##############################

def split_suffix(n):
	if n < 1 or n > 26:
		raise ValueError(f"Only 1–26 split genes supported; got {n}")
	return string.ascii_lowercase[n - 1]

##############################

def read_gff(path):
	records = []
	genes = {}
	transcripts = {}
	tx_by_gene = defaultdict(list)
	child_features_by_tx = defaultdict(list)

	print(f"Reading GFF from {path}", file=sys.stderr )
	with open(path) as fh:
		for line in fh:
			if line.startswith("#") or not line.strip():
				records.append(("raw", line))
				continue

			parts = line.rstrip("\n").split("\t")
			if len(parts) != 9:
				records.append(("raw", line))
				continue

			attrs = parse_attrs(parts[8])
			feat = {
				"parts": parts,
				"seqid": parts[0],
				"source": parts[1],
				"type": parts[2],
				"start": int(parts[3]),
				"end": int(parts[4]),
				"score": parts[5],
				"strand": parts[6],
				"phase": parts[7],
				"attrs": attrs,
				"line": line.rstrip("\n"),
			}

			records.append(("feature", feat))

			if feat["type"] == "gene" and "ID" in attrs:
				genes[attrs["ID"]] = feat

			elif feat["type"] in {"transcript", "mRNA"} and "ID" in attrs:
				tx_id = attrs["ID"]
				transcripts[tx_id] = feat
				for p in parent_ids(attrs):
					tx_by_gene[p].append(tx_id)

			else:
				for p in parent_ids(attrs):
					child_features_by_tx[p].append(feat)
	print("Found {} genes and {} transcripts".format( len(genes), len(transcripts) ), file=sys.stderr )

	return records, genes, transcripts, tx_by_gene, child_features_by_tx

##############################

def load_splits_table(path):
	splits = defaultdict(lambda: defaultdict(list))

	print(f"Reading false fusions table from {path}", file=sys.stderr )
	with open(path) as fh:
		header = next(fh).rstrip("\n").split("\t")
		idx_gene = header.index("old_gene_id")
		idx_cluster = header.index("cluster")
		idx_tx = header.index("transcript_id")

		for line in fh:
			if not line.strip():
				continue
			fields = line.rstrip("\n").split("\t")
			gene_id = fields[idx_gene]
			cluster = int(fields[idx_cluster])
			tx_id = fields[idx_tx]
			splits[gene_id][cluster].append(tx_id)

	final = {}
	for gene_id, clusters in splits.items():
		final[gene_id] = [
			clusters[i]
			for i in sorted(clusters)
		]
	print("Found {} fusion candidate genes".format( len(splits) ), file=sys.stderr )
	return final

##############################

def parse_false_exon_string(s):
	"""
	Format:
	  c1:4348-5664:-
	"""
	seqid, rest = s.split(":", 1)
	coords, strand = rest.rsplit(":", 1)
	start, end = coords.split("-")
	return seqid, int(start), int(end), strand

##############################

def load_motif_summary(path):
	"""
	Reads summary TSV from detect script.

	Required columns:
	  transcript_id
	  false_first_exon
	  remove_first_exon

	Only rows with remove_first_exon=yes are used for deletion.
	"""
	linecounter = 0
	false_exon_by_tx = {}
	print("Reading motif matches table from {}".format( path ), file=sys.stderr )
	with open(path) as fh:
		header = next(fh).rstrip("\n").split("\t")
		idx_tx = header.index("transcript_id")
		idx_fe = header.index("false_first_exon")
		idx_remove = header.index("remove_first_exon")

		for line in fh:
			linecounter += 1
			if not line.strip():
				continue

			fields = line.rstrip("\n").split("\t")

			if fields[idx_remove] != "yes":
				continue

			tx_id = fields[idx_tx]
			fe = fields[idx_fe]

			if fe:
				false_exon_by_tx[tx_id] = parse_false_exon_string(fe)
	print("Found {} motifs, and {} flagged for removal".format( linecounter, len(false_exon_by_tx) ), file=sys.stderr )
	return false_exon_by_tx

##############################

def build_rename_maps(splits):
	old_gene_to_new_genes = defaultdict(list)
	old_tx_to_new_gene = {}
	old_tx_to_new_tx = {}

	for old_gene_id, clusters in splits.items():
		for gene_i, tx_cluster in enumerate(clusters, start=1):
			new_gene_id = f"{old_gene_id}{split_suffix(gene_i)}"
			old_gene_to_new_genes[old_gene_id].append(new_gene_id)

			for tx_i, old_tx_id in enumerate(tx_cluster, start=1):
				new_tx_id = f"{new_gene_id}.i{tx_i}"
				old_tx_to_new_gene[old_tx_id] = new_gene_id
				old_tx_to_new_tx[old_tx_id] = new_tx_id

	return old_gene_to_new_genes, old_tx_to_new_gene, old_tx_to_new_tx

##############################

def is_false_first_exon(feat, old_tx_id, false_exon_by_tx):
	if feat["type"] != "exon":
		return False

	if old_tx_id not in false_exon_by_tx:
		return False

	seqid, start, end, strand = false_exon_by_tx[old_tx_id]

	return (
		feat["seqid"] == seqid
		and feat["start"] == start
		and feat["end"] == end
		and feat["strand"] == strand
	)

##############################

def get_new_tx_id(old_tx_id, old_tx_to_new_tx):
	return old_tx_to_new_tx.get(old_tx_id, old_tx_id)

##############################

def get_new_gene_id_for_tx(old_tx_id, transcripts, old_tx_to_new_gene):
	if old_tx_id in old_tx_to_new_gene:
		return old_tx_to_new_gene[old_tx_id]

	old_parent = parent_ids(transcripts[old_tx_id]["attrs"])[0]
	return old_parent

##############################

def compute_corrected_spans(
	transcripts,
	tx_by_gene,
	child_features_by_tx,
	false_exon_by_tx,
	old_tx_to_new_gene,
	old_tx_to_new_tx,
):
	tx_spans = {}
	gene_spans = {}

	old_tx_to_final_gene = {}

	for old_tx_id, tx in transcripts.items():
		if old_tx_id in old_tx_to_new_gene:
			final_gene = old_tx_to_new_gene[old_tx_id]
		else:
			parents = parent_ids(tx["attrs"])
			if not parents:
				continue
			final_gene = parents[0]

		old_tx_to_final_gene[old_tx_id] = final_gene

		kept_children = []

		for feat in child_features_by_tx.get(old_tx_id, []):
			if is_false_first_exon(feat, old_tx_id, false_exon_by_tx):
				continue
			kept_children.append(feat)

		exon_like = [
			f for f in kept_children
			if f["type"] in {"exon", "CDS", "five_prime_UTR", "three_prime_UTR", "UTR"}
		]

		if exon_like:
			start = min(f["start"] for f in exon_like)
			end = max(f["end"] for f in exon_like)
		else:
			# If removing the first exon leaves no child features,
			# keep the original transcript span rather than deleting transcript.
			start = tx["start"]
			end = tx["end"]

		final_tx = get_new_tx_id(old_tx_id, old_tx_to_new_tx)

		tx_spans[old_tx_id] = {
			"new_tx_id": final_tx,
			"gene_id": final_gene,
			"seqid": tx["seqid"],
			"start": start,
			"end": end,
			"strand": tx["strand"],
		}

		if final_gene not in gene_spans:
			gene_spans[final_gene] = {
				"seqid": tx["seqid"],
				"start": start,
				"end": end,
				"strand": tx["strand"],
			}
		else:
			gene_spans[final_gene]["start"] = min(gene_spans[final_gene]["start"], start)
			gene_spans[final_gene]["end"] = max(gene_spans[final_gene]["end"], end)

	return tx_spans, gene_spans, old_tx_to_final_gene

##############################

def write_repaired_gff(
	records,
	genes,
	transcripts,
	false_exon_by_tx,
	splits,
	old_gene_to_new_genes,
	old_tx_to_new_gene,
	old_tx_to_new_tx,
	tx_spans,
	gene_spans,
	output_gff,
):
	split_gene_ids = set(splits)
	written_split_genes = set()

	def write_split_gene(out, old_gene_id, new_gene_id):
		old_gene = genes[old_gene_id]
		span = gene_spans[new_gene_id]

		new_attrs = dict(old_gene["attrs"])
		new_attrs["ID"] = new_gene_id
		new_attrs["Name"] = new_gene_id
		new_attrs["Original_ID"] = old_gene_id

		parts = list(old_gene["parts"])
		parts[0] = span["seqid"]
		parts[3] = str(span["start"])
		parts[4] = str(span["end"])
		parts[6] = span["strand"]
		parts[8] = format_attrs(new_attrs)

		out.write("\t".join(parts) + "\n")
		written_split_genes.add(new_gene_id)

	with open(output_gff, "w") as out:
		for kind, obj in records:
			if kind == "raw":
				out.write(obj)
				continue

			feat = obj
			ftype = feat["type"]
			attrs = dict(feat["attrs"])

			if ftype == "gene":
				old_gene_id = attrs.get("ID")

				if old_gene_id in split_gene_ids:
					# Do not write all split gene records here.
					# Each split gene is written immediately before the
					# first transcript assigned to it.
					continue

				if old_gene_id in gene_spans:
					span = gene_spans[old_gene_id]
					attrs["Name"] = attrs.get("Name", old_gene_id)

					parts = list(feat["parts"])
					parts[3] = str(span["start"])
					parts[4] = str(span["end"])
					parts[8] = format_attrs(attrs)
					out.write("\t".join(parts) + "\n")
					continue

				out.write(feat["line"] + "\n")
				continue

			if ftype in {"transcript", "mRNA"}:
				old_tx_id = attrs.get("ID")

				if old_tx_id in tx_spans:
					span = tx_spans[old_tx_id]

					new_tx_id = span["new_tx_id"]
					new_gene_id = span["gene_id"]

					# If this transcript belongs to a split gene, write that
					# split gene immediately before the first transcript for it.
					if old_tx_id in old_tx_to_new_gene:
						old_gene_id = feat["attrs"].get("Parent", "")

						if (
							new_gene_id not in written_split_genes
							and old_gene_id in genes
						):
							write_split_gene(out, old_gene_id, new_gene_id)

					if old_tx_id in old_tx_to_new_tx:
						attrs["Original_ID"] = old_tx_id
						attrs["Original_Parent"] = feat["attrs"].get("Parent", "")

					attrs["ID"] = new_tx_id
					attrs["Parent"] = new_gene_id

					parts = list(feat["parts"])
					parts[3] = str(span["start"])
					parts[4] = str(span["end"])
					parts[8] = format_attrs(attrs)

					out.write("\t".join(parts) + "\n")
					continue

				out.write(feat["line"] + "\n")
				continue

			parents = parent_ids(attrs)

			if parents:
				# Drop motif-hit first exon from all affected transcripts.
				drop_feature = False
				for old_parent in parents:
					if is_false_first_exon(feat, old_parent, false_exon_by_tx):
						drop_feature = True
						break

				if drop_feature:
					continue

				new_parents = []
				changed = False

				for old_parent in parents:
					new_parent = get_new_tx_id(old_parent, old_tx_to_new_tx)
					new_parents.append(new_parent)
					if new_parent != old_parent:
						changed = True

				if changed:
					attrs["Parent"] = ",".join(new_parents)

					if "ID" in attrs:
						old_feature_id = attrs["ID"]
						attrs["Original_ID"] = old_feature_id

						for old_tx_id, new_tx_id in old_tx_to_new_tx.items():
							if old_tx_id in old_feature_id:
								attrs["ID"] = old_feature_id.replace(
									old_tx_id,
									new_tx_id,
									1,
								)
								break

					parts = list(feat["parts"])
					parts[8] = format_attrs(attrs)

					out.write("\t".join(parts) + "\n")
					continue

			out.write(feat["line"] + "\n")

##############################

def write_summary(path, false_exon_by_tx, old_tx_to_new_tx, old_tx_to_new_gene, tx_spans):
	with open(path, "w") as out:
		out.write(
			"old_tx_id\tnew_tx_id\tnew_gene_id\tremoved_first_exon\t"
			"new_tx_start\tnew_tx_end\n"
		)

		for old_tx_id in sorted(false_exon_by_tx):
			fe = false_exon_by_tx[old_tx_id]
			fe_string = f"{fe[0]}:{fe[1]}-{fe[2]}:{fe[3]}"

			span = tx_spans.get(old_tx_id, {})
			out.write(
				f"{old_tx_id}\t"
				f"{old_tx_to_new_tx.get(old_tx_id, old_tx_id)}\t"
				f"{old_tx_to_new_gene.get(old_tx_id, span.get('gene_id', ''))}\t"
				f"{fe_string}\t"
				f"{span.get('start', '')}\t"
				f"{span.get('end', '')}\n"
			)

##############################

def main():
	ap = argparse.ArgumentParser()
	ap.add_argument("--gff", required=True)
	ap.add_argument("--splits-table", required=True)
	ap.add_argument("--motif-summary", required=True)
	ap.add_argument("--out-gff", required=True)
	ap.add_argument("--repair-summary", required=True)
	args = ap.parse_args()

	records, genes, transcripts, tx_by_gene, child_features_by_tx = read_gff(args.gff)

	splits = load_splits_table(args.splits_table)
	false_exon_by_tx = load_motif_summary(args.motif_summary)

	old_gene_to_new_genes, old_tx_to_new_gene, old_tx_to_new_tx = build_rename_maps(splits)

	tx_spans, gene_spans, old_tx_to_final_gene = compute_corrected_spans(
		transcripts=transcripts,
		tx_by_gene=tx_by_gene,
		child_features_by_tx=child_features_by_tx,
		false_exon_by_tx=false_exon_by_tx,
		old_tx_to_new_gene=old_tx_to_new_gene,
		old_tx_to_new_tx=old_tx_to_new_tx,
	)

	write_repaired_gff(
		records=records,
		genes=genes,
		transcripts=transcripts,
		false_exon_by_tx=false_exon_by_tx,
		splits=splits,
		old_gene_to_new_genes=old_gene_to_new_genes,
		old_tx_to_new_gene=old_tx_to_new_gene,
		old_tx_to_new_tx=old_tx_to_new_tx,
		tx_spans=tx_spans,
		gene_spans=gene_spans,
		output_gff=args.out_gff,
	)

	write_summary(
		path=args.repair_summary,
		false_exon_by_tx=false_exon_by_tx,
		old_tx_to_new_tx=old_tx_to_new_tx,
		old_tx_to_new_gene=old_tx_to_new_gene,
		tx_spans=tx_spans,
	)

	print(f"Wrote repaired GFF: {args.out_gff}")
	print(f"Wrote repair summary: {args.repair_summary}")


if __name__ == "__main__":
	main()

example_results = """
Reading GFF from Hcv1av93.gff
Found 14591 genes and 20076 transcripts
Reading false fusions table from Hcv1av93_proposed_clusters.tsv
Found 7 fusion candidate genes
Reading motif matches table from Hcv1av93_splits_summary.tsv
Found 598 motifs, and 499 flagged for removal
Wrote repaired GFF: Hcv1av93.leader_removed.gff
Wrote repair summary: Hcv1av93.leader_removed.repair_summary.tsv

Reading GFF from Hcv1av93_v10.gff
Found 14619 genes and 19660 transcripts
Reading false fusions table from Hcv1av93_v10_proposed_clusters.tsv
Found 7 fusion candidate genes
Reading motif matches table from Hcv1av93_v10_splits_summary.tsv
Found 601 motifs, and 501 flagged for removal
Wrote repaired GFF: Hcv1av93_v10.leader_removed.gff
Wrote repair summary: Hcv1av93_v10.leader_removed.repair_summary.tsv
"""
