#!/usr/bin/env python
#
# get_read_skip_from_bam.py

'''
get_read_skip_from_bam.py  last modified 2026-06-04
    from a SAM or BAM file
    extract length of skipped bases at beginning or end of a long read

~/samtools-1.9/samtools view UCSC_Hcal_v1_B1_LR.sorted.bam | get_read_skip_from_bam.py - > UCSC_Hcal_v1_B1_LR.sorted.read_skip.tab

    creates a 4-column table, of line number, starting -S, ending -S, and sequence of -S

cut -f 4 UCSC_Hcal_v1_B1_LR.sorted.read_skip.tab | sort | uniq -c | sort -nr

1845104 0
 354257 CCCTCAAAGTTTGAAAAGTTGTGATGAAATTTGTTTAATTAAAC
 321762 GGGAGTTTCAAACTTTTCAACACTACTTTAAACAAATTAATTTG

'''

import sys
import argparse
import re
import time
from collections import defaultdict

def cigar_to_list(cigarstring):
	"""from a cigar string, return a list of each number letter pair
    e.g. cigar_to_list("1S347M97N85M2388N79M107N175M110N176M111N227M118N100M89N102M124N272M110N1393M1D271M1I4M39S")
    should return
    ['1S', '347M', '97N', '85M', '2388N', '79M', '107N', '175M', '110N', '176M', '111N', '227M', '118N', '100M', '89N', '102M', '124N', '272M', '110N', '1393M', '1D', '271M', '1I', '4M', '39S']
    """
	cigar_list = [] # build growing list with each RE
	for rematch in re.finditer(r"\d+",cigarstring):
		rematch_start = rematch.start() # start of numbers \d+
		rematch_end = rematch.end() + 1 # should always be 1 letter after one or more numbers
		number_letter_pair = cigarstring[rematch_start:rematch_end]
		cigar_list.append(number_letter_pair)
	return cigar_list


def use_skipped_sequence(motif, polyN_cutoff):
	'''check if motif is polyA or polyT, and return boolean'''
	motif_len = len(motif)
	T_percent = motif.count("T") * 1.0 / motif_len
	A_percent = motif.count("A") * 1.0 / motif_len
	if T_percent >= polyN_cutoff or A_percent >= polyN_cutoff:
		return False
	else:
		return True


def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('input_file', type = argparse.FileType('r'), default = '-', help="SAM file, - for stdin from BAM")
	parser.add_argument('-n', '--target-min', type=int, default=30, help="print out leader sequences with S value of at least N [30]", metavar="N")
	parser.add_argument('-N', '--target-max', type=int, default=50, help="print out leader sequences with S value of at most N [50]", metavar="N")
	parser.add_argument('-r', '--repeat-percent', type=float, default=0.9, help="cutoff to ignore polyA or polyT [0.9]", metavar="0.N")
	args = parser.parse_args(argv)

	linecounter = 0
	cigar_counter = 0

	match_target_count = 0
	polyA_count = 0
	no_cigar_count = 0

	hardclip_count = defaultdict(int)
	verbose = True

	sys.stderr.write("# Reading {}, tracking S of {}-{}bp  {}\n".format(args.input_file.name, args.target_min, args.target_max, time.asctime() ) )
	for line in args.input_file:
		linecounter += 1
		lsplits = line.split("\t")
		# QNAME FLAG RNAME POS  MAPQ CIGAR RNEXT PNEXT  TLEN SEQ QUAL extra
		if len(lsplits) < 12:
			continue
		scaffold = lsplits[2]
		cigar_string = lsplits[5]
		cigar_list = cigar_to_list(cigar_string)
		if len(cigar_list)==0:
			no_cigar_count += 1
			continue

		if cigar_list[0][-1]=="S": # for 5prime
			five_pr_skip = int(cigar_list[0][0:-1])
		else:
			five_pr_skip = 0
		if cigar_list[-1][-1]=="S": # for 3prime
			three_pr_skip = int(cigar_list[-1][0:-1])
		else:
			three_pr_skip = 0
		cigar_counter += 1
		if cigar_list[0][-1]=="H" or cigar_list[-1][-1]=="H":
			hardclip_count[scaffold] += 1

		leader_seq = "0"

		seq_string = lsplits[9]
		# get sequence at 5prime if within length
		if args.target_min <= five_pr_skip <= args.target_max:
			match_target_count += 1
			skipped_seq = seq_string[0:five_pr_skip]
			if use_skipped_sequence(skipped_seq, args.repeat_percent):
				leader_seq = skipped_seq
			else:
				polyA_count += 1
		# get sequence at 3prime
		if args.target_min <= three_pr_skip <= args.target_max:
			match_target_count += 1
			seq_string_len = len(seq_string)
			skipped_seq = seq_string[ (seq_string_len - three_pr_skip):]
			if use_skipped_sequence(skipped_seq, args.repeat_percent):
				leader_seq = skipped_seq
			else:
				polyA_count += 1

		sys.stdout.write("{}\t{}\t{}\t{}\n".format( linecounter, five_pr_skip, three_pr_skip, leader_seq ) )

	sys.stderr.write("# Counted {} lines for {} reads  {}\n".format(linecounter, cigar_counter, time.asctime() ) )
	if no_cigar_count:
		sys.stderr.write("# {} had no CIGAR string\n".format( no_cigar_count ) )
	if match_target_count:
		sys.stderr.write("# {} reads had either end matched to target leader length\n".format(match_target_count) )
	if polyA_count:
		sys.stderr.write("# {} were likely polyA tails, and were ignored\n".format( polyA_count ) )

	# number should sum to the same number of hardclips in the bam file
	#
	# samtools view SRR10403849_vs_UCSC_Hcal_v1.bam | awk '$6 ~ /H/ {print $2}' | sort | uniq -c

	if hardclip_count and verbose is True:
		sys.stderr.write("# {} hard clips were found ({:.2f}%)\n".format( sum(hardclip_count.values()), 100.8*sum(hardclip_count.values())/cigar_counter ) )
		for k,v in hardclip_count.items():
			print("{}\t{}".format( k,v ), file=sys.stderr)

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
