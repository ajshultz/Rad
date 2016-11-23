#! /usr/bin/env python

import re, sys, os, itertools, sets, getopt        # Load standard modules I often use

'''
This script will take the output from an MS simulation and output a file with the distribution of Fishers exact test p-values for allele frequency differences between the populations.
'''

def main(argv):
	try:
		opts,args = getopt.getopt(argv,'hio',)
	except getopt.GetOptError:
		print "ms_output_pvals.py -h help -i <input file (ms results file)> -o <output file."
		sys.exit(2)
		
	input = ''
	output = 'ms_pvals_out.csv'
	
	for opt, arg in opts:
		if opt == "-h":
			print "ms_output_pvals.py -h help -i <input file (ms results file)> -o <output file."
			sys.exit(2)
		elif opt == "-i":
			input = arg
		elif opt == "-o":
			output = arg
			

	class MSRecord(object):
	
		def __init__(self,
	
	
	
	
	
	ms_out = open(input,"r")
	
	for line in ms_out:
		if line[0] == "m"