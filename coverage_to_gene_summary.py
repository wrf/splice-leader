#!/usr/bin/env python3
# v1.0 2026-05-11

"""coverage_to_gene_summary.py  v1.0 2026-05-11

  used in pipeline with gff_exons_to_unique_bed.py

bedtools coverage -s -split -a reference_exons.bed -b set1.bam > set1.coverage.tsv
bedtools coverage -s -split -a reference_exons.bed -b set2.bam > set2.coverage.tsv

    such as:

gff_exons_to_unique_bed.py -g GCF_026151205.1_MBARI_Bmic_1.0_genomic.gff -o GCF_026151205.1_MBARI_Bmic_1.0_genomic.exons.bed

bedtools coverage -s -split -a GCF_026151205.1_MBARI_Bmic_1.0_genomic.exons.bed -b BolinopsisDeep-Iso-V4656_vs_MBARI_Bmic_1.0.bam > BolinopsisDeep-Iso-V4656_vs_MBARI_Bmic_1.0.coverage.tsv
bedtools coverage -s -split -a GCF_026151205.1_MBARI_Bmic_1.0_genomic.exons.bed -b BolinopsisDeep-Iso-V4656_vs_MBARI_Bmic_1.0.w_leader.bam > BolinopsisDeep-Iso-V4656_vs_MBARI_Bmic_1.0.w_leader.coverage.tsv

~/git/splice-leader/coverage_to_gene_summary.py -q BolinopsisDeep-Iso-V4656_vs_MBARI_Bmic_1.0.coverage.tsv -s BolinopsisDeep-Iso-V4656_vs_MBARI_Bmic_1.0.w_leader.coverage.tsv -o BolinopsisDeep-Iso-V4656_vs_MBARI_Bmic_1.0.w_leader.coverage.tab

"""

import sys
import argparse
from collections import defaultdict

############################################################

def parse_gene_id(name_field):
	"""
	Extract gene_id from BED name column.

	Assumes:
	  - genbank IDs for exons as gene=XYZ
	"""
	if "gene=" in name_field:
		for part in name_field.split(";"):
			if part.startswith("gene="):
				return part.split("=")[1]
	return name_field

############################################################

def read_coverage_file(path):
	"""
	Reads bedtools coverage output.

	Returns list of dicts:
	{
		chrom, start, end, name, strand,
		count, covered_bases, length, fraction
	}
	"""
	records = []
	print("# Reading coverage table from {}".format(path), file=sys.stderr)
	with open(path) as f:
		for line in f:
			if not line.strip():
				continue

			fields = line.rstrip("\n").split("\t")

			chrom = fields[0]
			start = int(fields[1])
			end = int(fields[2])
			name = fields[3]
			strand = fields[5]

			# bedtools coverage appended fields (last 4)
			count = int(fields[-4])
			covered_bases = int(fields[-3])
			length = int(fields[-2])
			fraction = float(fields[-1])

			records.append({
				"chrom": chrom,
				"start": start,
				"end": end,
				"name": name,
				"strand": strand,
				"count": count,
				"covered_bases": covered_bases,
				"length": length,
				"fraction": fraction
			})
	print("# Found {} records from {}".format( len(records), path), file=sys.stderr)
	return records

############################################################

def exon_is_covered(record, min_count=1, min_fraction=0.1):
	"""
	Define coverage rule.

	You can tune this.
	"""
	return (record["count"] >= min_count) and (record["fraction"] >= min_fraction)

############################################################

def main():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument("-q", "--cov1", required=True, help="Coverage TSV for set1")
	parser.add_argument("-s", "--cov2", required=True, help="Coverage TSV for set2")
	parser.add_argument("-o", "--output", default="-")
	parser.add_argument("--min-count", type=int, default=1)
	parser.add_argument("--min-fraction", type=float, default=0.1)

	args = parser.parse_args()

	cov1 = read_coverage_file(args.cov1)
	cov2 = read_coverage_file(args.cov2)

	if len(cov1) != len(cov2):
		raise ValueError("Coverage files must have same number of rows and same ordering")

	print("# Compiling coverage, using min count of {}, minimum fraction of {}".format(args.min_count, args.min_fraction), file=sys.stderr)
	# gene_id -> list of exon tuples
	gene_to_exons = defaultdict(list)

	for r1, r2 in zip(cov1, cov2):
		gene_id = parse_gene_id(r1["name"])

		exon = {
			"chrom": r1["chrom"],
			"start": r1["start"],
			"end": r1["end"],
			"strand": r1["strand"],
			"cov1": r1,
			"cov2": r2
		}

		gene_to_exons[gene_id].append(exon)

	out = sys.stdout if args.output == "-" else open(args.output, "w")

	header = [
		"gene_id", "scaffold", "start_pos", "end_pos", "span",
		"n_reference_exons",
		"n_matched_exons1", "fraction_matched1", "bitstring1",
		"n_matched_exons2", "fraction_matched2", "status_flag"
	]
	out.write("\t".join(header) + "\n")

	nomatch_total = 0
	nls_count = 0

	for gene_id in sorted(gene_to_exons):
		exons = gene_to_exons[gene_id]

		# sort exons
		exons = sorted(exons, key=lambda x: (x["chrom"], x["start"], x["end"]))

		starts = [e["start"] for e in exons]
		ends = [e["end"] for e in exons]

		scaffold = exons[0]["chrom"]
		gene_start = min(starts)
		gene_end = max(ends)
		span = gene_end - gene_start

		bits1 = []
		bits2 = []

		matched1 = 0
		matched2 = 0

		for exon in exons:
			covered1 = exon_is_covered(exon["cov1"], args.min_count, args.min_fraction)
			covered2 = exon_is_covered(exon["cov2"], args.min_count, args.min_fraction)

			bits1.append("1" if covered1 else "0")
			bits2.append("1" if covered2 else "0")

			if covered1:
				matched1 += 1
			if covered2:
				matched2 += 1

		total = len(exons)

		frac1 = matched1 / total if total else 0.0
		frac2 = matched2 / total if total else 0.0

		if matched1 == 0 and matched2 == 0:
			nomatch_total += 1
		status_flag = "-"
		if matched1 > 0 and matched2 == 0:
			status_flag = "nls"
			nls_count += 1

		out.write(
			f"{gene_id}\t{scaffold}\t{gene_start}\t{gene_end}\t{span}\t"
			f"{total}\t{matched1}\t{frac1:.2f}\t{''.join(bits1)}\t"
			f"{matched2}\t{frac2:.2f}\t{status_flag}\n"
		)

	if out is not sys.stdout:
		out.close()
	print("# {} genes no matches in either query".format(nomatch_total), file=sys.stderr)
	print("# {} genes had matches in query but not secondary".format(nls_count), file=sys.stderr)

############################################################

if __name__ == "__main__":
	main()
