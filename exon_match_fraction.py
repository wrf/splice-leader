#!/usr/bin/env python3
# v1 created 2026-03-31

"""
exon_match_fraction.py  v1.1

For each reference gene, compute the fraction of its exon intervals
that are matched by at least one exon in the query GFF, using
interval-based matching rather than exact coordinate equality.
"""

import sys
import argparse
from collections import defaultdict
from bisect import bisect_right

##############################

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
			else:
				attrs[part.strip()] = ""
		return attrs

	parts = [x.strip() for x in attr_text.split(";") if x.strip()]
	for part in parts:
		if " " in part:
			key, value = part.split(" ", 1)
			attrs[key.strip()] = value.strip().strip('"')
		else:
			attrs[part.strip()] = ""
	return attrs

##############################

def get_feature_id(attrs):
	for key in ("ID", "gene_id", "transcript_id", "Name"):
		if key in attrs and attrs[key]:
			return attrs[key]
	return None

##############################

def get_parent_ids(attrs):
	for key in ("Parent",):
		if key in attrs and attrs[key]:
			return [x.strip() for x in attrs[key].split(",") if x.strip()]
	if "transcript_id" in attrs and attrs["transcript_id"]:
		return [attrs["transcript_id"]]
	if "gene_id" in attrs and attrs["gene_id"]:
		return [attrs["gene_id"]]
	return []

##############################

def normalize_strand(strand):
	if strand in ("+", "-"):
		return strand
	return "."

##############################

def parse_reference_gff(reference_gff, ignore_strand=False):
	gene_ids = set()
	transcript_to_gene = {}
	gene_to_exons = defaultdict(set)

	print("# Reading reference GFF from {}".format(reference_gff), file=sys.stderr)

	with open(reference_gff, "r", encoding="utf-8") as fh:
		for line in fh:
			if not line.strip() or line.startswith("#"):
				continue

			fields = line.rstrip("\n").split("\t")
			if len(fields) != 9:
				continue

			scaffold, source, feature_type, start, end, score, strand, phase, attr_text = fields
			attrs = parse_attributes(attr_text)
			strand = normalize_strand(strand)

			if feature_type == "gene":
				fid = get_feature_id(attrs)
				if fid:
					gene_ids.add(fid)

			elif feature_type in ("mRNA", "transcript"):
				tx_id = get_feature_id(attrs)
				parents = get_parent_ids(attrs)
				if tx_id and parents:
					transcript_to_gene[tx_id] = parents[0]

	print("# Found {} unique gene IDs from {}".format(len(gene_ids), reference_gff), file=sys.stderr)

	with open(reference_gff, "r", encoding="utf-8") as fh:
		for line in fh:
			if not line.strip() or line.startswith("#"):
				continue

			fields = line.rstrip("\n").split("\t")
			if len(fields) != 9:
				continue

			scaffold, source, feature_type, start, end, score, strand, phase, attr_text = fields
			if feature_type != "exon":
				continue

			attrs = parse_attributes(attr_text)
			strand = normalize_strand(strand)
			start = int(start)
			end = int(end)

			if ignore_strand:
				exon_key = (scaffold, start, end)
			else:
				exon_key = (scaffold, start, end, strand)

			parents = get_parent_ids(attrs)
			if not parents:
				continue

			assigned_gene_ids = set()

			for parent in parents:
				if parent in transcript_to_gene:
					assigned_gene_ids.add(transcript_to_gene[parent])
				elif parent in gene_ids:
					assigned_gene_ids.add(parent)
				else:
					if "gene_id" in attrs and attrs["gene_id"]:
						assigned_gene_ids.add(attrs["gene_id"])

			for gid in assigned_gene_ids:
				gene_to_exons[gid].add(exon_key)

	print("# Found {} genes with exons from {}".format(len(gene_to_exons), reference_gff), file=sys.stderr)
	return gene_to_exons

##############################

def parse_query_exons(query_gff, ignore_strand=False):
	index = defaultdict(list)

	print("# Reading query GFF from {}".format(query_gff), file=sys.stderr)

	with open(query_gff, "r", encoding="utf-8") as fh:
		for line in fh:
			if not line.strip() or line.startswith("#"):
				continue

			fields = line.rstrip("\n").split("\t")
			if len(fields) != 9:
				continue

			scaffold, source, feature_type, start, end, score, strand, phase, attr_text = fields
			if feature_type != "exon":
				continue

			start = int(start)
			end = int(end)
			strand = normalize_strand(strand)

			if ignore_strand:
				key = scaffold
			else:
				key = (scaffold, strand)

			index[key].append((start, end))

	final_index = {}
	total_intervals = 0

	for key, intervals in index.items():
		intervals = sorted(set(intervals))
		starts = [x[0] for x in intervals]
		final_index[key] = {
			"intervals": intervals,
			"starts": starts
		}
		total_intervals += len(intervals)

	print("# Indexed {} unique query exons from {}".format(total_intervals, query_gff), file=sys.stderr)
	return final_index

##############################

def get_gene_bounds(exon_set, ignore_strand=False):
	if not exon_set:
		return None

	if ignore_strand:
		seqids = {e[0] for e in exon_set}
		if len(seqids) != 1:
			raise ValueError("Exons span multiple seqids")

		seqid = next(iter(seqids))
		starts = [e[1] for e in exon_set]
		ends = [e[2] for e in exon_set]

		return (seqid, min(starts), max(ends))

	else:
		seqids = {e[0] for e in exon_set}
		strands = {e[3] for e in exon_set}

		if len(seqids) != 1:
			raise ValueError("Exons span multiple seqids")
		if len(strands) != 1:
			raise ValueError("Exons span multiple strands")

		seqid = next(iter(seqids))
		strand = next(iter(strands))

		starts = [e[1] for e in exon_set]
		ends = [e[2] for e in exon_set]

		return (seqid, min(starts), max(ends), strand)

##############################

def compute_all_gene_bounds(gene_to_exons, ignore_strand=False):
	gene_to_bounds = {}
	for gene_id, exon_set in gene_to_exons.items():
		bounds = get_gene_bounds(exon_set, ignore_strand=ignore_strand)
		gene_to_bounds[gene_id] = bounds
	print("# Computed gene bounds for {} genes".format(len(gene_to_bounds)), file=sys.stderr)
	return gene_to_bounds

##############################

def exon_length(start, end):
	return end - start + 1

##############################

def overlap_bp(a_start, a_end, b_start, b_end):
	ov_start = max(a_start, b_start)
	ov_end = min(a_end, b_end)
	if ov_start > ov_end:
		return 0
	return ov_end - ov_start + 1

##############################

def ref_coverage_fraction(ref_start, ref_end, q_start, q_end):
	ov = overlap_bp(ref_start, ref_end, q_start, q_end)
	ref_len = ref_end - ref_start + 1
	return ov / ref_len if ref_len > 0 else 0.0

##############################

def get_candidate_intervals(exon, query_index, ignore_strand=False):
	if ignore_strand:
		scaffold, start, end = exon
		key = scaffold
	else:
		scaffold, start, end, strand = exon
		key = (scaffold, strand)

	if key not in query_index:
		return []

	intervals = query_index[key]["intervals"]
	starts = query_index[key]["starts"]

	right = bisect_right(starts, end)

	candidates = []
	for q_start, q_end in intervals[:right]:
		if q_end >= start:
			candidates.append((q_start, q_end))

	return candidates

##############################

def exon_has_interval_match(exon, query_index, ignore_strand=False,
							min_ref_coverage=0.8,
							boundary_tolerance=1):
	if ignore_strand:
		scaffold, ref_start, ref_end = exon
	else:
		scaffold, ref_start, ref_end, strand = exon

	candidates = get_candidate_intervals(exon, query_index, ignore_strand=ignore_strand)

	for q_start, q_end in candidates:
		ref_cov = ref_coverage_fraction(ref_start, ref_end, q_start, q_end)
		start_shift = abs(ref_start - q_start)
		end_shift = abs(ref_end - q_end)

		if ref_cov >= min_ref_coverage:
			return True

		if start_shift <= boundary_tolerance and end_shift <= boundary_tolerance:
			return True

	return False

##############################

def exon_match_bitstring(ref_exons, query_index, ignore_strand=False,
						 min_ref_coverage=0.8,
						 boundary_tolerance=1):
	ordered = sorted(ref_exons, key=lambda x: (x[0], x[1], x[2], x[3] if len(x) > 3 else "."))

	bits = []
	bit_sum = 0
	for exon in ordered:
		if exon_has_interval_match(
			exon, 
			query_index,
			ignore_strand=ignore_strand,
			min_ref_coverage=min_ref_coverage,
			boundary_tolerance=boundary_tolerance
		):
			bits.append("1")
			bit_sum += 1
		else:
			bits.append("0")

	return "".join(bits), bit_sum

############################################################

def main():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument("-r", "--reference", required=True, help="Reference GFF/GTF")
	parser.add_argument("-q", "--query", required=True, help="Query GFF/GTF")
	parser.add_argument("-s", "--secondary", help="Secondary query GFF/GTF")
	parser.add_argument("-o", "--output", default="-", help="Output TSV file (default: stdout)")
	parser.add_argument("--ignore-strand", action="store_true", help="Match exons by scaffold/start/end only, ignoring strand")
	parser.add_argument("--min-ref-coverage", type=float, default=0.8, help="Minimum fraction of the reference exon covered by a query exon [default: 0.8]")
	parser.add_argument("--boundary-tolerance", type=int, default=1, help="Allow exon boundaries to differ by this many bp [default: 1]")
	args = parser.parse_args()

	ref_gene_to_exons = parse_reference_gff(args.reference, ignore_strand=args.ignore_strand)
	gene_bounds = compute_all_gene_bounds(ref_gene_to_exons, ignore_strand=args.ignore_strand)
	query_exons = parse_query_exons(args.query, ignore_strand=args.ignore_strand)
	secondary_exons = parse_query_exons(args.secondary, ignore_strand=args.ignore_strand) if args.secondary else {}

	out_fh = sys.stdout if args.output == "-" else open(args.output, "w", encoding="utf-8")
	print("# Output will be written to {}".format(out_fh), file=sys.stderr)

	header = [
		"gene_id", "scaffold", "start_pos", "end_pos", "span",
		"n_reference_exons",
		"n_matched_exons1", "fraction_matched1", "bitstring1",
		"n_matched_exons2", "fraction_matched2", "status_flag"
	]
	out_fh.write("\t".join(header) + "\n")

	ref_exon_total = 0
	query1_match_total = 0
	query2_match_total = 0
	nomatch_total = 0
	q_to_s_nomatch_total = 0

	for gene_id in sorted(ref_gene_to_exons):
		ref_exons = ref_gene_to_exons[gene_id]
		total_exons = len(ref_exons)
		ref_exon_total += total_exons
		bounds = gene_bounds[gene_id]
		scaffold = bounds[0]
		gene_span = bounds[2] - bounds[1] + 1
		no_leader_flag = "-"

		bitstring1, matched1 = exon_match_bitstring(
			ref_exons,
			query_exons,
			ignore_strand=args.ignore_strand,
			min_ref_coverage=args.min_ref_coverage,
			boundary_tolerance=args.boundary_tolerance
		)

		if args.secondary:
			bitstring2, matched2 = exon_match_bitstring(
				ref_exons,
				secondary_exons,
				ignore_strand=args.ignore_strand,
				min_ref_coverage=args.min_ref_coverage,
				boundary_tolerance=args.boundary_tolerance
			)
		else:
			bitstring2 = "0" * total_exons
			matched2 = 0

		query1_match_total += matched1
		query2_match_total += matched2

		fraction1 = matched1 / total_exons if total_exons else 0.0
		fraction2 = matched2 / total_exons if total_exons else 0.0

		if matched1 == 0 and matched2 == 0:
			nomatch_total += 1
		if matched1 > 0 and matched2 == 0:
			q_to_s_nomatch_total += 1
			no_leader_flag = "nls"

		out_fh.write(
			f"{gene_id}\t{scaffold}\t{bounds[1]}\t{bounds[2]}\t{gene_span}\t"
			f"{total_exons}\t{matched1}\t{fraction1:.2f}\t{bitstring1}\t"
			f"{matched2}\t{fraction2:.2f}\t{no_leader_flag}\n"
		)

	if out_fh is not sys.stdout:
		out_fh.close()

	print("# Found {} query matches out of {} reference exons".format(query1_match_total, ref_exon_total), file=sys.stderr)
	print("# Found {} secondary matches out of {} reference exons".format(query2_match_total, ref_exon_total), file=sys.stderr)
	print("# {} genes no matches in either query".format(nomatch_total), file=sys.stderr)
	print("# {} genes had matches in query but not secondary".format(q_to_s_nomatch_total), file=sys.stderr)

##############################

if __name__ == "__main__":
	main()
