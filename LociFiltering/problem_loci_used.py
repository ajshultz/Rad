#! /usr/bin/env python

import re, sys, os, itertools, sets, getopt        # Load standard modules I often use

"""
This script will take the loci from a vcf file, and identify any that have been identified as a problem with catalog_read_pairs.py output (either both 1 and 2 reads, or mapping to multiple loci).  t\Two files will be produced in the output director: Blacklist.txt which contains a list with the name of each problem locus on a single line, and ProblemStats.txt, which is a tab delimited file containing stats on how many loci were problematic within individuals, across individuals, in either categories, and the number of  loci not in either of those categories.

Note that the input directory needs to contain the "MultipleLociWithinInds.txt" and "MultipleLociAcrossInds.txt" files created from catalog_read_pair.py.
"""

def main(argv):
	try:
		opts,args = getopt.getopt(argv,'hv:w:a:b:',)
	except getopt.GetOptError:
		print "problem_loci_used.py -v <Path to vcf file> -o <path to output directory (default .)> -i <path to input directory with problem files created in catalog_read_pair.py>"
		sys.exit(2)
			
	vcffile = ''
	outputdir = '.'
	inputdir = '.'
	
	for opt, arg in opts:
		if opt == "-h":
			print "problem_loci_used.py -v <Path to vcf file> -o <path to output directory (default .)> -i <path to input directory with problem files created in catalog_read_pair.py>"
			sys.exit(2)
		elif opt == "-v":
			vcffile = arg
		elif opt == "-i":
			inputdir = arg
		elif opt == "-o":
			outputdir = arg

	blacklist_file = "%s/Blacklist.txt"%outputdir
	stats_file = "%s/ProblemStats.txt"%outputdir
	within_ind_file = '%s/MultipleLociWithinInds.csv'%inputdir
	across_ind_file ='%s/MultipleLociAcrossInds.csv'%inputdir	
	

	within_ind = open(within_ind_file,"r")
	across_ind = open(across_ind_file,"r")
	vcf = open(vcffile,"r")
	stats = open(stats_file,"w")
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

	stats.write("Bad loci within ind\t%d\n"%(len(set(poploci) & set(withinlist))))
	stats.write("Bad loci across ind\t%d\n"%(len(set(poploci) & set(acrosslist))))
	stats.write("Bad loci in either\t%d\n"%(len(bad_used)))
	stats.write("Loci with single matches\t%d"%(len(poploci) - len(bad_used)))

	blacklist = open(blacklist_file,"w")
	for locus in bothlist:
		blacklist.write("%s\n"%locus)
	
	within_ind.close()
	across_ind.close()
	vcf.close()
	blacklist.close()
	
if __name__ == "__main__":
		main(sys.argv[1:])