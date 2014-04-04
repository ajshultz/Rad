#! /usr/bin/env python

import re, sys, os, itertools, sets, getopt        # Load standard modules I often use
import Bio
from Bio import Blast
from Bio.Blast import NCBIXML

"""
This script will take in a blast xml file and see how many hits each locus had.  It will output a table of the number of hits, the number of loci with each hit number, and the percentage of loci with that number of hits.  Also, if a blacklist file is added (produced from problem_loci_used.py, the program will calculate the number of reads with more than 1 hit detected by the pairing method.
If the flag -n is added, a new blacklist file will be created.
"""

def main(argv):
	try:
		opts,args = getopt.getopt(argv,'hb:o:p:n:',)
	except getopt.GetOptError:
		print "BlastResultsNumMatches.py -b <path to blast xml file> -o <path to output file (default ./numhits.csv)> -p <path to blacklist file> -n <path to new blacklist file (default = none)>"
		sys.exit(2)
			
	blastfile=''
	outputfile = './numhits.csv'
	blacklistfile = 0
	newblacklist = 0

	for opt, arg in opts:
		if opt == "-h":
			print "BlastResultsNumMatches.py -b <path to blast xml file> -o <path to output file (default ./numhits.csv)> -p <path to blacklist file> -n <path to new blacklist file (default = none)>"
			sys.exit(2)
		elif opt == "-b":
			blastfile = arg
		elif opt == "-o":
			outputfile = arg
		elif opt == "-p":
			blacklistfile = arg
		elif opt == "-n":
			newblacklist = arg



	#First, open the Blast output file, and make a dictionary of how many hits occurred for each locus, and a dictionary with the snpid (from the consensus sequence position) as the key and the chromosome and base pair positions (start and stop) as the value.

	blast = open(blastfile,"r")
	blast_records = Bio.Blast.NCBIXML.parse(blast)
	blast_dict = {}
	for record in blast_records:
		blast_dict[record.query]=record


	#Make a function for and extract the number of alignments for each locus - into a new dictionary.

	def NumBlastAlignments(blast_dict,tag_id):
		record = blast_dict[str(tag_id)]
		numaligns = 0
		for alignment in record.alignments:
			for hsp in alignment.hsps:
				numaligns += 1
		return numaligns

	loci =  blast_dict.keys()
	numhitsdict = {}
	hitstatsdict = {}

	for locus in loci:
		numhitsdict[locus] = NumBlastAlignments(blast_dict,locus)
		if numhitsdict[locus] in hitstatsdict:
			hitstatsdict[numhitsdict[locus]] += 1
		else:
			hitstatsdict[numhitsdict[locus]] = 1
	
	#Create a dictionary with the number of hits as the key and a list of the loci with that number as the value
	hitslocidict = {}
	for hits in hitstatsdict.keys():
		hitslocidict[hits] = [locus for locus,numaligns in numhitsdict.iteritems() if numaligns == hits]
	
	if newblacklist:
		newblacklistout = open(newblacklist,"w")
		for hits in hitslocidict:
			if hits > 1:
				for locus in hitslocidict[hits]:
					newblacklistout.write("%s\n"%(locus))
		newblacklistout.close()
		

	hitstats = sorted(zip(hitstatsdict.keys(),hitstatsdict.values()))

	if blacklistfile:
		blacklist = open(blacklistfile,"r")
		output = open(outputfile,"w")
		output.write("Number of hits,Number of loci,Percent of loci,Num loci in blacklist,Percent of category in blacklist,Num left after blacklist\n")
		problem = []
		for line in blacklist:
			line = line.strip()
			problem.append(line)
		for hit in hitstats:
			percent = float(hit[1])/float(len(blast_dict.keys()))
			inproblem = []
			for locus in hitslocidict[hit[0]]:
				if locus in problem:
					inproblem.append(locus)
				else:
					pass
			numinproblem = len(inproblem)
			percinproblem = numinproblem/float(hit[1])
			afterblacklist = hit[1]-numinproblem
			output.write("%d,%d,%f,%d,%f,%d\n"%(hit[0],hit[1],percent,numinproblem,percinproblem,afterblacklist))				
		blacklist.close()
	else:		
		output = open(outputfile,"w")
		output.write("Number of hits,Number of loci,Percent of loci\n")
		for hit in hitstats:
			percent = float(hit[1])/float(len(blast_dict.keys()))
			output.write("%d,%d,%f\n"%(hit[0],hit[1],percent))


	blast.close()
	output.close()
	

if __name__ == "__main__":
		main(sys.argv[1:])