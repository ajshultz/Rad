#! /usr/bin/env python

import re, sys, os, itertools, sets, getopt        # Load standard modules I often use
import pickle

"""
This script will take the loci from a vcf file, and identify any that have been identified as a problem (either both 1 and 2 reads, or mapping to multiple loci).
"""

def main(argv):
	try:
		opts,args = getopt.getopt(argv,'hv:w:a:b:',)
	except getopt.GetOptError:
		print "problem_loci_used.py -v <Path to vcf file> -w <Path to problem loci within individuals> -a <Path to problem loci across individuals> -b <Directory for blacklist results>"
		sys.exit(2)
			
	vcffile = ''
	within_ind_file = ''
	across_ind_file =''
	blacklist_dir = '.'

	for opt, arg in opts:
		if opt == "-h":
			print "catalog_read_pair.py -l <mysql host (default localhost)> -u <mysql user (default root)> -p <mysql password> -s <Stacks database name> -o <path to output directory (default .)> -n <number of samples> -m <minimum cut off depth for a stack to be included in catalogs (default 1)>"
			sys.exit(2)
		elif opt == "-v":
			vcffile = arg
		elif opt == "-w":
			within_ind_file = arg
		elif opt == "-a":
			across_ind_file = arg
		elif opt == "-b":
			blacklist_dir = arg

	blacklist_file = "%s/Blacklist.txt"%blacklist_dir


	within_ind = open(within_ind_file,"r")
	across_ind = open(across_ind_file,"r")
	vcf = open(vcffile,"r")
	withinlist = []

	for line in within_ind:
		line = line.strip()
		withinlist.append(line)
	
	acrosslist = []

	for line in across_ind:
		newline = line.split(",")
		acrosslist.append(newline[0])
	
	bothlist = set(withinlist) | set(acrosslist)
	print len(bothlist)

	poploci = []

	for line in vcf:
		if line[0] != "#":
			newline = line.split("\t")
			if newline[2] not in poploci:
				poploci.append(newline[2])
			else:
				pass
		else:
			pass

	bad_used = set(bothlist) & set(poploci)

	print "Bad loci within ind", len(set(poploci) & set(withinlist))
	print "Bad loci across ind", len(set(poploci) & set(acrosslist))
	print "Bad loci in either", len(bad_used)

	print len(poploci) - len(bad_used)

	blacklist = open(blacklist_file,"w")
	for locus in bothlist:
		blacklist.write("%s\n"%locus)
	
	within_ind.close()
	across_ind.close()
	vcf.close()
	blacklist.close()
	
if __name__ == "__main__":
		main(sys.argv[1:])