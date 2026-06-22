#!/usr/bin/env python
# leader_sequence_uniqc_to_fasta.py v1.0 created 2025-12-27
# v1.1 add reverse complement 2026-05-04

"""leader_sequence_uniqc_to_fasta.py  v1.1 2026-05-04

Convert lines like:
  31770 GAGTTTTAATACTTTCAACACTACTATATAACAAATAATTTGAGG

(from: uniq -c | sort -nr)
into FASTA:

>prefix|count
sequence

>Bolinopsis|31770
GAGTTTTAATACTTTCAACACTACTATATAACAAATAATTTGAGG

NOTE: names of sequences will be unique, but the sequences may not be

Usage examples:
gzip -dc Bolinopsis_Isoseq.fasta.gz | grep -v ">" | cut -c 1-45 | sort | uniq -c | sort -nr | head -n 100 | leader_sequence_uniqc_to_fasta.py --prefix Bolinopsis > Bolinopsis_Isoseq.leader_sequences.fasta

if the uniq -c output was already captured, this can be accepted by piping:
cat leader_seq.counts | ~/git/splice-leader/leader_sequence_uniqc_to_fasta.py -m 100 --prefix Bmicroptera
"""

import sys
import argparse

_RC_TABLE = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')

def parse_line(line: str):
	"""
	Parse a single 'uniq -c' line.
	Returns (count:int, item:str) or (None, None) if not parseable.
	"""
	s = line.strip("\n")
	if not s.strip():
		return None, None

	# Split only on first whitespace after the count
	parts = s.lstrip().split(None, 1)
	if len(parts) < 2:
		return None, None

	count_str, item = parts[0], parts[1]
	try:
		count = int(count_str)
	except ValueError:
		return None, None

	item = item.rstrip()
	if not item:
		return None, None

	return count, item

def write_fasta(out, header: str, seq: str):
	out.write(f">{header}\n")
	out.write(seq + "\n")

def revcomp(seq: str) -> str:
    return seq.translate(_RC_TABLE)[::-1]

def main():
	ap = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	ap.add_argument('-m', "--minimum-count", type=int, default=1, help="Skip entries whose count is less than N [1].")
	ap.add_argument("--prefix", default="",
					help="Optional prefix for FASTA headers (e.g., 'Species_name').")
	args = ap.parse_args()

	needs_revcomp_CTENO = ["CAA","CTC",  "TTA","TCA","TTT","TTG","TGT","TAA","TAG",  "ATT","AAT","ATA","AAA","ATG","AGT"]

	n_lines = 0
	n_written = 0
	n_written_sum = 0
	revcomp_count = 0
	min_remove_count = 0 # n sequences removed
	min_remove_sum = 0 # total occurrences of the removed N sequences
	short_error_count = 0
	for line in sys.stdin:
		count, item = parse_line(line)
		n_lines += 1
		if count is None: # if error in parse_line() step
			continue
		if count < args.minimum_count: # below -m count
			min_remove_count += 1
			min_remove_sum += count
			continue
		if len(item) < 3: # seq was mostly empty - sometimes A T C or G alone
			short_error_count += 1
			continue
		header = f"{args.prefix}|{n_lines}|{count}"

		if item[0:3] in needs_revcomp_CTENO: # reverse complement some sequences
			revcomp_count += 1
			item = revcomp(item)
			header += "|rc"

		write_fasta(sys.stdout, header, item)
		n_written += 1
		n_written_sum += count

	print(f"Read {n_lines} lines, wrote {n_written} FASTA records, with total of {n_written_sum}.", file=sys.stderr)
	if revcomp_count:
		print(f"Found {revcomp_count} sequences likely should be reverse complement.", file=sys.stderr)
	if min_remove_count:
		min_remove_pct = 100.0*min_remove_sum/n_written_sum
		print(f"Skipped {min_remove_count} sequences with count less than {args.minimum_count}, for {min_remove_sum} total ({min_remove_pct:.2f}%).", file=sys.stderr)
	if short_error_count:
		print(f"Ignored {short_error_count} sequences under length of 3bp.", file=sys.stderr)

if __name__ == "__main__":
	main()

