#!/usr/bin/env python3
# v1.0 2026-05-11

"""gff_exons_to_unique_bed.py  v1.0  last modified 2026-05-11

    used in pipeline with:

~/git/splice-leader/gff_exons_to_unique_bed.py -g GCF_026151205.1_MBARI_Bmic_1.0_genomic.gff -o GCF_026151205.1_MBARI_Bmic_1.0_genomic.exons.bed
bedtools coverage -s -split -a GCF_026151205.1_MBARI_Bmic_1.0_genomic.exons.bed -b BolinopsisDeep-Iso-V4656_vs_MBARI_Bmic_1.0.w_leader.bam > BolinopsisDeep-Iso-V4656_vs_MBARI_Bmic_1.0.w_leader.coverage.tsv
~/git/splice-leader/coverage_to_gene_summary.py -q BolinopsisDeep-Iso-V4656_vs_MBARI_Bmic_1.0.coverage.tsv -s BolinopsisDeep-Iso-V4656_vs_MBARI_Bmic_1.0.w_leader.coverage.tsv -o BolinopsisDeep-Iso-V4656_vs_MBARI_Bmic_1.0.w_leader.coverage.tab

"""

import sys
import argparse

def parse_attributes(attr_text):
	attrs = {}
	attr_text = attr_text.strip()
	if not attr_text:
		return attrs

	if "=" in attr_text:
		parts = [x.strip() for x in attr_text.split(";") if x.strip()]
		for part in parts:
			if "=" in part:
				key, value = part.split("=", 1)
				attrs[key.strip()] = value.strip().strip('"')
		return attrs

	parts = [x.strip() for x in attr_text.split(";") if x.strip()]
	for part in parts:
		if " " in part:
			key, value = part.split(" ", 1)
			attrs[key.strip()] = value.strip().strip('"')
	return attrs


def get_feature_id(attrs):
	for key in ("ID", "gene_id", "transcript_id", "Name"):
		if key in attrs and attrs[key]:
			return attrs[key]
	return None


def get_parent_ids(attrs):
	if "Parent" in attrs and attrs["Parent"]:
		return [x.strip() for x in attrs["Parent"].split(",") if x.strip()]
	if "transcript_id" in attrs and attrs["transcript_id"]:
		return [attrs["transcript_id"]]
	if "gene_id" in attrs and attrs["gene_id"]:
		return [attrs["gene_id"]]
	return []


def normalize_strand(strand):
	return strand if strand in ("+", "-") else "."


def main():
	parser = argparse.ArgumentParser(
		description="Convert GFF/GTF exon features to BED, collapsing identical exons within each gene."
	)
	parser.add_argument("-g", "--gff", required=True, help="Reference GFF/GTF")
	parser.add_argument("-o", "--output", default="-", help="Output BED file")
	args = parser.parse_args()

	gene_ids = set()
	transcript_to_gene = {}

	# pass 1: collect genes and transcript->gene mapping
	print("# Reading GFF from {}".format(args.gff), file=sys.stderr)
	with open(args.gff, "r", encoding="utf-8") as fh:
		for line in fh:
			if not line.strip() or line.startswith("#"):
				continue
			fields = line.rstrip("\n").split("\t")
			if len(fields) != 9:
				continue

			chrom, source, feature_type, start, end, score, strand, phase, attr_text = fields
			attrs = parse_attributes(attr_text)

			if feature_type == "gene":
				fid = get_feature_id(attrs)
				if fid:
					gene_ids.add(fid)

			elif feature_type in ("mRNA", "transcript"):
				tx_id = get_feature_id(attrs)
				parents = get_parent_ids(attrs)
				if tx_id and parents:
					transcript_to_gene[tx_id] = parents[0]
	print("# Found {} gene IDs and {} transcripts from {}".format( len(gene_ids), len(transcript_to_gene), args.gff), file=sys.stderr)

	# pass 2: collect unique exons by gene+coords+strand
	unique_exons = set()
	print("# Reading unique exons from {}".format(args.gff), file=sys.stderr)
	with open(args.gff, "r", encoding="utf-8") as fh:
		for line in fh:
			if not line.strip() or line.startswith("#"):
				continue
			fields = line.rstrip("\n").split("\t")
			if len(fields) != 9:
				continue

			chrom, source, feature_type, start, end, score, strand, phase, attr_text = fields
			if feature_type != "exon":
				continue

			start = int(start)
			end = int(end)
			strand = normalize_strand(strand)
			attrs = parse_attributes(attr_text)
			parents = get_parent_ids(attrs)

			assigned_gene_ids = set()

			for parent in parents:
				if parent in transcript_to_gene:
					assigned_gene_ids.add(transcript_to_gene[parent])
				elif parent in gene_ids:
					assigned_gene_ids.add(parent)

			if not assigned_gene_ids and "gene_id" in attrs and attrs["gene_id"]:
				assigned_gene_ids.add(attrs["gene_id"])

			for gene_id in assigned_gene_ids:
				# BED is 0-based half-open
				bed_start = start - 1
				bed_end = end
				unique_exons.add((chrom, bed_start, bed_end, gene_id, strand))
	print("# Found {} unique exons".format( len(unique_exons) ), file=sys.stderr)
	sorted_exons = sorted(unique_exons, key=lambda x: (x[0], x[3], x[1], x[2], x[4]))

	out = sys.stdout if args.output == "-" else open(args.output, "w", encoding="utf-8")
	for chrom, bed_start, bed_end, gene_id, strand in sorted_exons:
		out.write(f"{chrom}\t{bed_start}\t{bed_end}\t{gene_id}\t.\t{strand}\n")

	if out is not sys.stdout:
		out.close()


if __name__ == "__main__":
	main()
