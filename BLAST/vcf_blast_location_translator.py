#! /usr/bin/env python

import re, sys, os, itertools, sets, getopt        # Load standard modules I often use
import MySQLdb   # Gives the ability to interact with the mysql database.
import pickle
import vcf
import Bio
from Bio import Blast
from Bio.Blast import NCBIXML

"""
This script will take in a vcf file output from stacks, a csv hit table output from blasts for the corresponding tags, and a mysql database.  It will access the snp position for the snp ids present in the vcf file from the plink map file, and will translate them to the chromosomal positions from the blast output.  Any tags that blasted to more than one location on the zebra finch genome or did not have any hits will be dropped from the vcf file. 

The input: -v <path to vcf file to be translated> -b <path to blast xml file> -o <path to output file (e.g. ./output.vcf)> -m <path to plink map file> -t <path to desired tag id output file (e.g. ./output.tag_id.txt)> -s <if 1, only a single snp per locus will be written - if 0 (default) all loci will be translated>

If singleSnp = 1, only one snp per locus will be written.

A file with all tags numbers in the VCF file will also be written.
"""


def main(argv):
	try:
		opts,args = getopt.getopt(argv,'hv:b:o:m:t:s:',)
	except getopt.GetOptError:
		print "vcf_blast_location_translator.py -v <path to vcf file to be translated> -b <path to blast xml file> -o <path to output file (e.g. ./output.vcf)> -m <path to plink map file> -t <path to desired tag id output file (e.g. ./output.tag_id.txt)> -s <if 1, only a single snp per locus will be written - if 0 (default) all loci will be translated>"
		sys.exit(2)
			
	vcffile = ''
	mapfile= ''
	output = ''
	blastfile = ''
	singleSnp = 0
	tagfile = ''


	for opt, arg in opts:
		if opt == "-h":
			print "vcf_blast_location_translator.py -v <path to vcf file to be translated> -b <path to blast xml file> -o <path to output file> -m <path to plink map file> -t <path to desired tag id output file> -s <if 1, only a single snp per locus will be written - if 0 (default) all loci will be translated>"
			sys.exit(2)
		elif opt == "-v":
			vcffile = arg
		elif opt == "-b":
			blastfile = arg
		elif opt == "-m":
			mapfile = arg
		elif opt == "-o":
			output = arg
		elif opt == "-t":
			tagfile = arg
		elif opt == "-s":
			singleSnp = arg


	#chrtransfile = "/Users/allisonshultz/Dropbox/PythonScripts/RAD/BLAST/ChromsomeNameConvBlast.csv"

	#First, open the Blast output file, and make a dictionary of how many hits occurred for each locus, and a dictionary with the snpid (from the consensus sequence position) as the key and the chromosome and base pair positions (start and stop) as the value.

	blast = open(blastfile,"r")
	blast_records = Bio.Blast.NCBIXML.parse(blast)
	blast_dict = {}
	for record in blast_records:
		blast_dict[record.query]=record

	#Make a chromosome translation dictionary - not longer necessary for the moment, extracting this directly out of blast
	# chrtrans = open(chrtransfile,"r")
	# chrtransdict = {}
	# for line in chrtrans:
	# 	line = line.strip().split(",")
	# 	chrtransdict[line[1]] = line[0]
	# 
	# chrtrans.close()



	#Open the input and output VCF files, write all commented data to new file unaltered

	vcfin = open(vcffile,"r")
	vcfout = open(output,"w")
	map = open(mapfile,"r")


	#This function takes in a blast output dictionary that is parsed form NCBIXML, a tag_id and position along the Stacks input fragment.  It then returns the corresponding chromsome and snp position in the Zebra Finch genome as a list.  Note that in some cases BLAST query was shorter than the Stacks fragment.  In the case where the SNP was not in the region that mapped to the ZF genome, the function will return 0.
	def NewSnpPos(blast_dict,tag_id,fragment_pos):
		record = blast_dict[str(tag_id)]
		pos = int(fragment_pos)-1#We have to subtract 1 from pos to allow for python subsetting which starts at 0
		for alignment in record.alignments:
			for hsp in alignment.hsps:
				if pos > hsp.query_start -1 + len(hsp.query) - hsp.query.count("-"):
					return 0,0
				else:
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
		for description in record.descriptions:
			chr = str(description.title.strip().split(",")[0].split()[-1])
		return chr,sbjct_snppos

	def NumBlastAlignments(blast_dict,tag_id):
		record = blast_dict[str(tag_id)]
		numaligns = 0
		for alignment in record.alignments:
			for hsp in alignment.hsps:
				numaligns += 1
		return numaligns

	#print NewSnpPos(blast_dict,524,1)
	#print NumBlastAlignments(blast_dict,9665)


	#Create a dictionary of the snpid as the key and the tag number _ snp postion on the fragment as the value.  Note that I am finding that the vcf and plink snpids are one off from one another, and am correcting.  I believe this was fixed in a more recent version of stacks, but appears to be off again... so if this script throws an error with a different vcf file that may be the problem.
	maploci = {}
	for line in map:
		if line[0] == "#":
			pass
		else:
			info = line.strip().split("\t")
			maploci[info[3]] = info[1]

	uniquetag = []

	for line in vcfin:
		if line[0] == "#":
			vcfout.write(line)
		else:
			vcfline = line.strip().split("\t")
			snpid = vcfline[1]
			tag_pos = maploci[str(int(snpid) - 1)]#This is where the correction is.
			#tag_pos = maploci[snpid]
			tag,pos = tag_pos.split("_")
			if singleSnp == 1:
				if tag not in uniquetag:
					uniquetag.append(tag)
					if NumBlastAlignments(blast_dict,tag) == 1:
						chr,newpos = NewSnpPos(blast_dict,tag,pos)
						if newpos != 0:#In case they mapped beyond the mapped blast read
							if chr != "scaffold":#So that we only get the results that actually mapped to the chromosomes.
								vcfline[0] = str(chr)
								vcfline[1] = str(newpos)
								newvcfline = "\t".join(vcfline)
								vcfout.write("%s\n"%(newvcfline))
			else:
				if NumBlastAlignments(blast_dict,tag) == 1:
					chr,newpos = NewSnpPos(blast_dict,tag,pos)
					if newpos != 0:#In case they mapped beyond the mapped blast read
						if chr != "scaffold":#So that we only get the results that actually mapped to the chromosomes.
							vcfline[0] = str(chr)
							vcfline[1] = str(newpos)
							vcfline[2] = "%s_%s"%(tag,pos)
							print vcfline[2]
							newvcfline = "\t".join(vcfline)
							vcfout.write("%s\n"%(newvcfline))

	tags = open(tagfile,"w")
	for i in range(0,len(uniquetag)):
		tags.write("%s\n"%uniquetag[i])

	tags.close()	
	vcfin.close()
	vcfout.close()	
	map.close()
	blast.close()
	
	
if __name__ == "__main__":
		main(sys.argv[1:])