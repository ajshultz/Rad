#! /usr/bin/env python

import re, sys, os, itertools, sets, getopt        # Load standard modules I often use
import Bio
from Bio import Blast
from Bio.Blast import NCBIWWW

"""
This script takes in a fasta file, blasts it against the ZF genome with a given expect score, and returns an xml file.
"""

def main(argv):
	try:
		opts,args = getopt.getopt(argv,'hf:o:e:',)
	except getopt.GetOptError:
		print "BlastZFGenome.py -f <path to input fasta file> -o <path to output directory and filename (default ./blast.xml)> -e <expect value (default 10e-20)>" 
		sys.exit(2)
			
	output = './blast.xml'
	fastafile = ''
	expect = 10e-20

	for opt, arg in opts:
		if opt == "-h":
			print "BlastZFGenome.py -f <path to input fasta file> -o <path to output directory and filename (default ./blast.xml)> -e <expect value (default 10e-20)>"
			sys.exit(2)
		elif opt == "-o":
			output = arg
		elif opt == "-f":
			fastafile = arg
		elif opt == "-e":
			expect = float(arg)
				
	fasta = open(fastafile).read()
	result_handle = NCBIWWW.qblast("blastn","GPIPE/59729/101/ref_top_level",fasta,expect=expect)

	save_file = open(output, "w")
	save_file.write(result_handle.read())
	save_file.close()
	result_handle.close()

if __name__ == "__main__":
		main(sys.argv[1:])