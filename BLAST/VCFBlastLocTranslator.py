#! /usr/bin/env python

import re, sys, os, itertools, sets, getopt        # Load standard modules I often use
import MySQLdb   # Gives the ability to interact with the mysql database.
import pickle
import vcf

"""
This script will take in a vcf file output from stacks, a csv hit table output from blasts for the corresponding tags, and a mysql database.  It will access the snp position for the snp ids present in the vcf file from the plink map file, and will translate them to the chromosomal positions from the blast output.  Any tags that blasted to more than one location on the zebra finch genome or did not have any hits will be dropped from the vcf file. Note that if a file with chromosome name is provided, the ensemble chromosome names will be translated as well (chrname,ensembl name) for each line
"""

#User input variables
vcffile = "/Users/allisonshultz/Dropbox/HFRad-Tags/HFPaired_Reduced/_PhyloPop_r.5p14_allsitesoutput/batch_1.vcf"
mapfile= "/Users/allisonshultz/Dropbox/HFRad-Tags/HFPaired_Reduced/_PhyloPop_r.5p14_allsitesoutput/batch_1.plink.map"
output = "../TestResults/PhyloPop_r.5p14_translated.vcf"
blastfile = "/Users/allisonshultz/Dropbox/HFRad-Tags/PhyloPopr.5p14BLAST/HA4EU7HH014-Alignment-HitTable.csv"
chrtransfile = "/Users/allisonshultz/Dropbox/PythonScripts/RAD/BLAST/ChromsomeNameConvBlast.csv"

#First, open the Blast output file, and make a dictionary of how many hits occurred for each locus, and a dictionary with the snpid (from the consensus sequence position) as the key and the chromosome and base pair positions (start and stop) as the value.

blast = open(blastfile,"r")

locihits = {}
blastloc = {}
for line in blast:
	res = line.strip().split(",")
	if res[0] not in locihits:
		locihits[res[0]] = 1
		blastloc[res[0]] = [res[1],res[8],res[9]]
	else:
		locihits[res[0]] += 1

blast.close()

#Make a chromosome translation dictionary

chrtrans = open(chrtransfile,"r")

chrtransdict = {}
for line in chrtrans:
	line = line.strip().split(",")
	chrtransdict[line[1]] = line[0]


chrtrans.close()



#Open the input and output VCF files, write all commented data to new file unaltered

vcfin = open(vcffile,"r")
vcfout = open(output,"w")
map = open(mapfile,"r")

#Create a dictionary of the snpid as the key and the tag number _ snp postion on the fragment as the value.  Note that I am finding that the vcf and plink snpids are one off from one another, and am correcting.  I believe this was fixed in a more recent version of stacks, so if this script throws and error with a different vcf file that may be the problem.

maploci = {}
for line in map:
	info = line.strip().split("\t")
	maploci[info[3]] = info[1]

for line in vcfin:
	if line[0] == "#":
		vcfout.write(line)
	else:
		vcfline = line.strip().split("\t")
		snpid = vcfline[1]
		tag_pos = maploci[str(int(snpid) - 1)]#This is where the correction is.
		tag,pos = tag_pos.split("_")
		if tag in blastloc:
			numhits = locihits[tag]
			if numhits == 1:		
				zfloc = blastloc[tag]
				zfchr = chrtransdict[zfloc]
				zfstart = 
			else:
				pass
		else:
			pass

	
vcfin.close()
vcfout.close()	
map.close()