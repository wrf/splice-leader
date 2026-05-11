#!/usr/bin/env python
# created 2025-11-26

"""
transfer_gff_annotation.py ML2.3.21_annotations.gff ML2.3.vs.Mlei_v2.gff > ML2.3.vs.Mlei_v2.w_annot.gff
"""

import sys
import time

#ML0001	ML2.3.21	gene	17958	18941	.	-	.	ID=ML0001.g4;Name=ML0001.g4;description=homeobox transcription factor HD77b;user_note=annotated from Ryan2010
#ML0001	ML2.3.21	transcript	17958	18941	.	-	.	ID=ML0001.g4.i1;Parent=ML0001.g4;source_ID=Mlei_GFAT01_c358047_g2_i2
#ML0001	ML2.3.21	exon	17958	18941	.	-	.	Parent=ML0001.g4.i1

#CM108477.1	pinfish	transcript	17726	18220	.	-	.	ID=ML0572.g1.i1;Name=ML0572.g1.i1;
#CM108477.1	pinfish	exon	17726	18220	.	-	.	Parent=ML0572.g1.i1

# get(feature,'strand')==1? '#005824' : '#66c2a4'


if len(sys.argv) < 2:
	sys.exit( __doc__ )
else:
	line_counter = 0
	annot_dict = {} # key is gene ID, value is Description and all other tags
	sys.stderr.write("# Parsing annotation GFF {}  {}\n".format(sys.argv[1], time.asctime() ) )
	for line in open(sys.argv[1],'r'):
		if line[0]=="#":
			continue
		line_counter += 1
		lsplits = line.split('\t')
		feature = lsplits[2]
		if feature not in ["gene","transcript"]:
			continue
		attributes = lsplits[8].strip()
		attrd = dict([(field.strip().split("=",1)) for field in attributes.split(";") if field.count("=")])
		extra_annot_string = ""
		gene_id = attrd["ID"]
		for item in attrd.items():
			if item[0] not in ["ID","Parent","Name"]:
				extra_annot_string += "{}={};".format(*item)
		if extra_annot_string:
			annot_dict[gene_id] = extra_annot_string
	sys.stderr.write("# Found {} GFF lines, {} descriptions  {}\n".format(line_counter, len(annot_dict), time.asctime() ) )

	line_counter = 0
	sys.stderr.write("# Parsing new GFF {}  {}\n".format(sys.argv[2], time.asctime() ) )
	for line in open(sys.argv[2],'r'):
		if line[0]=="#":
			sys.stdout.write(line)
			continue
		line_counter += 1
		lsplits = line.split('\t')
		feature = lsplits[2]
		if feature == "exon": # other are transcript
			sys.stdout.write(line)
			continue
		attributes = lsplits[8].strip()
		attrd = dict([(field.strip().split("=",1)) for field in attributes.split(";") if field.count("=")])
		gene_id = attrd.get("ID")
		parent_id = attrd.get("ID").rsplit('.',1)[0]
		if annot_dict.get(gene_id,False):
			attributes += annot_dict.get(gene_id)
		if annot_dict.get(parent_id,False):
			attributes += annot_dict.get(parent_id)
		lsplits[8] = attributes
		print( "\t".join(lsplits) , file=sys.stdout )
	sys.stderr.write("# Found {} GFF lines  {}\n".format(line_counter, time.asctime() ) )


