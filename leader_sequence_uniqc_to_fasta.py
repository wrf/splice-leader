#!/usr/bin/env python
# leader_sequence_uniqc_to_fasta.py v1.0 created 2025-12-27

"""leader_sequence_uniqc_to_fasta.py  v1.0 2025-12-27

Convert lines like:
  <count><whitespace><sequence...>

(from: uniq -c | sort -nr)
into FASTA:

>prefix|count
sequence

Usage examples:
gzip -dc Bolinopsis_Isoseq.fasta.gz | grep -v ">" | cut -c 1-45 | sort | uniq -c | sort -nr | head -n 40 | leader_sequence_uniqc_to_fasta.py --prefix Bolinopsis > Bolinopsis_Isoseq.leader_sequences.fasta
"""

import sys
import argparse

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

def main():
	ap = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	ap.add_argument("--skip-zero", action="store_true",
					help="Skip entries whose count is 0.")
	ap.add_argument("--prefix", default="",
					help="Optional prefix for FASTA headers (e.g., 'Species_name').")
	args = ap.parse_args()

	n_written = 0
	for line in sys.stdin:
		count, item = parse_line(line)
		if count is None:
			continue
		if args.skip_zero and count == 0:
			continue

		header = f"{args.prefix}|{count}"
		write_fasta(sys.stdout, header, item)
		n_written += 1

	print(f"Wrote {n_written} FASTA records.", file=sys.stderr)

if __name__ == "__main__":
	main()

