#!/usr/bin/env python3
# created by WRF 2025-12-02
# v1.1 2026-05-05 gzip output

"""retain_transcripts_w_leader_motif.py v1.1  last modified 2026-05-05

  Filter a FASTA file by motifs at the start of each sequence.

Usage:
retain_transcripts_w_leader_motif.py -i input.fasta -m motifs.fasta -o output.fasta

  Motifs are contained in a fasta file, and are searched for an exact match at
  the start of each input sequence.
  For instance
>Hcal_354257_rc
GGGAGTTTCAAACTTTTCAACACTACTTTAAACAAATTAATTTG
"""

import sys
import time
import argparse
import gzip
from collections import defaultdict
from Bio import SeqIO

def main():
	if not len(sys.argv[1:]):
		sys.argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-i','--input', help="fasta format file, can be .gz")
	parser.add_argument('-m','--motifs', help="fasta file of motif sequences - motifs are expected to be short, 30-50bp")
	parser.add_argument('-o',"--output", help="Output FASTA file (filtered)")
	args = parser.parse_args()

	sys.stderr.write( "# Reading motifs from {}  {}\n".format( args.motifs, time.asctime() ) )
	motifs = SeqIO.to_dict(SeqIO.parse(args.motifs,"fasta"))
	if motifs:
		sys.stderr.write( "# Counted {} motifs from {}\n".format( len(motifs), args.motifs ) )
	else: # meaning none found
		raise SystemExit("Error: MOTIFS list is empty. Please check -m.")

	seqcount = 0
	retained_count = defaultdict(int)

	if args.output.endswith(".gz"): # determine gzip output, can be much slower
		open_output = gzip.open(args.output, "wt")  # text mode
		sys.stderr.write( "# Writing output to {} as gzip\n".format( args.output ) )
	else:
		open_output = open(args.output, "w")
		sys.stderr.write( "# Writing output to {}\n".format( args.output ) )

	with open_output as fout:
		# Parse records with SeqIO
		sys.stderr.write( "# Reading sequences from {}  {}\n".format( args.input, time.asctime() ) )
		if args.input.rsplit('.',1)[-1]=="gz":
			input_handler = gzip.open(args.input,'rt')
		else:
			input_handler = args.input
		for seq_record in SeqIO.parse(input_handler, "fasta"):
			seqcount += 1
			seq_str = str(seq_record.seq)
			seq_u = seq_str.upper()
			for m in motifs.values():
				m_up = str(m.seq).upper()
				if seq_u.startswith(m_up):
					retained_count[m.id] += 1
					SeqIO.write(seq_record, fout, "fasta")
					break # stop searching if a motif matches, since they SHOULD be unique
	sys.stderr.write( "# Read sequences from {}, wrote {} ({:.2f}%)  {}\n".format( seqcount, sum(retained_count.values()), sum(retained_count.values())*100.0/seqcount, time.asctime() ) )
	for k,v in retained_count.items():
		print("{}\t{}\t{}".format( k,str(motifs.get(k).seq),v ), file=sys.stdout)

if __name__ == "__main__":
	main()

