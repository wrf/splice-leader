#!/usr/bin/env python3
# created by WRF 2025-12-02
# v1.2 2026-05-25 added optional motif trimming mode
# v1.1 2026-05-05 gzip output

"""retain_transcripts_w_leader_motif.py v1.2  last modified 2026-05-25

  Filter a FASTA file by motifs at the start of each sequence.

Usage:
retain_transcripts_w_leader_motif.py -i input.fasta -m motifs.fasta -o output.fasta
retain_transcripts_w_leader_motif.py -i input.fasta -m motifs.fasta -o output.fasta --trim-motif

  Motifs are contained in a fasta file, and are searched for an exact match at
  the start of each input sequence.
  For instance
>Hcal_354257_rc
GGGAGTTTCAAACTTTTCAACACTACTTTAAACAAATTAATTTG

  By default, matching input sequences are written unchanged.

  With --trim-motif, matching input sequences are retained, but the leading
  motif sequence is removed before writing. If multiple motifs match the start
  of a sequence, the longest matching motif is used for trimming.
"""

import sys
import time
import argparse
import gzip
from collections import defaultdict
from Bio import SeqIO


def open_maybe_gzip(filename, mode="rt"):
	"""Open plain text or gzip-compressed files by filename suffix."""
	if filename.endswith(".gz"):
		return gzip.open(filename, mode)
	return open(filename, mode)


def main():
	if not len(sys.argv[1:]):
		sys.argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-i','--input', help="fasta format file, can be .gz")
	parser.add_argument('-m','--motifs', help="fasta file of motif sequences - motifs are expected to be short, 30-50bp")
	parser.add_argument('-o',"--output", help="Output FASTA file (filtered)")
	parser.add_argument('-t','--trim-motif', action='store_true', help="Retain matching sequences, but trim the matched leading motif before output")
	args = parser.parse_args()

	sys.stderr.write( "# Reading motifs from {}  {}\n".format( args.motifs, time.asctime() ) )
	with open_maybe_gzip(args.motifs, "rt") as motif_handle:
		motifs = SeqIO.to_dict(SeqIO.parse(motif_handle,"fasta"))
	if motifs:
		sys.stderr.write( "# Counted {} motifs from {}\n".format( len(motifs), args.motifs ) )
	else: # meaning none found
		raise SystemExit("Error: MOTIFS list is empty. Please check -m.")

	# Store motifs once as uppercase strings to avoid repeatedly converting them
	# for every input sequence.
	motif_items = []
	for m in motifs.values():
		motif_seq = str(m.seq).upper()
		if not motif_seq:
			print("Error: motif {} has an empty sequence - skipping".format(m.id), file=sys.stderr)
		motif_items.append((m.id, motif_seq))

	seqcount = 0
	retained_count = defaultdict(int)

	if args.output.endswith(".gz"): # determine gzip output, can be much slower
		open_output = gzip.open(args.output, "wt")  # text mode
		sys.stderr.write( "# Writing output to {} as gzip\n".format( args.output ) )
	else:
		open_output = open(args.output, "w")
		sys.stderr.write( "# Writing output to {}\n".format( args.output ) )

	if args.trim_motif:
		sys.stderr.write( "# Trimming matched leader motifs before output: {}\n".format( args.trim_motif ) )

	with open_output as fout:
		# Parse records with SeqIO
		sys.stderr.write( "# Reading sequences from {}  {}\n".format( args.input, time.asctime() ) )
		with open_maybe_gzip(args.input, 'rt') as input_handler:
			for seq_record in SeqIO.parse(input_handler, "fasta"):
				seqcount += 1
				seq_str = str(seq_record.seq)
				seq_u = seq_str.upper()

				matching_motifs = []
				for motif_id, motif_seq in motif_items:
					if seq_u.startswith(motif_seq):
						matching_motifs.append((motif_id, motif_seq))

				if not matching_motifs:
					continue

				# Original behavior effectively used the first matching motif. In trimming
				# mode, use the longest matching motif so nested/partial motifs trim to
				# the most specific leader sequence found.
				if args.trim_motif:
					matched_id, matched_seq = max(matching_motifs, key=lambda x: len(x[1]))
					trim_len = len(matched_seq)
					seq_record = seq_record[trim_len:]
				else:
					matched_id, matched_seq = matching_motifs[0]

				retained_count[matched_id] += 1
				SeqIO.write(seq_record, fout, "fasta")

	written = sum(retained_count.values())
	pct = written*100.0/seqcount if seqcount else 0.0
	sys.stderr.write( "# Read sequences from {}, wrote {} ({:.2f}%)  {}\n".format( seqcount, written, pct, time.asctime() ) )
	for k,v in retained_count.items():
		print("{}\t{}\t{}".format( k,str(motifs.get(k).seq),v ), file=sys.stdout)

if __name__ == "__main__":
	main()
