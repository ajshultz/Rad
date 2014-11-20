#! /usr/bin/env python

import Bio
from Bio import Blast
from Bio.Blast import NCBIXML

blast = open("/Users/allisonshultz/Dropbox/HFRad-Tags/PhyloPopr.5p14BLAST/HA4EU7HH014-Alignment_onlyfirst5.xml","r")

blast_records = Bio.Blast.NCBIXML.parse(blast)


pos = 134

pos = pos-1

record = next(blast_records)
#for record in blast_records:
for alignment in record.alignments:
	for hsp in alignment.hsps:
		start = hsp.query_start
		diff_from_1 = 1-start
		pos = pos+diff_from_1
		count = 0
		for nuc in hsp.query:
			if count <= pos:
				if nuc == "-":
					pos += 1
				else:
					pass
			else:
				pass
			count += 1
		#Now we have a position number that takes into account any indels and a start different than 1 in the query sequence and 
		numind = hsp.sbjct.count("-",0,pos)
		sbjct_pos = pos - numind
		sbjct_snppos = hsp.sbjct_start+sbjct_pos
		print sbjct_snppos





blast.close()