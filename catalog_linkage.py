#! /usr/bin/env python

import re, sys, os, itertools, sets, getopt        # Load standard modules I often use
import MySQLdb   # Gives the ability to interact with the mysql database.
import pickle
import vcf

"""
This script will obtain a list of locus catalog tag_id numbers from an inputted vcf file.
"""

#User input variables
vcffile = "/Users/allisonshultz/Dropbox/HFRad-Tags/HFPaired_Reduced/_PhyloPop_r.5p14/batch_1.vcf"
mysqlhost = "localhost"
mysqluser = "root"
mysqlpasswd = "hofi"
stacksdb = "HFdenovo_PairTestingReducedm3N3_radtags"
numsamples=78
outputdir = "TestResults"



vcf = open(vcffile,"r")

#Get list of loci from vcf file
lociused = []

for line in vcf:
	if line[0] != "#":
		newline = line.split("\t")
		if newline[2] not in lociused:
			lociused.append(newline[2])
		else:
			pass
	else:
		pass


#Create database connection
MyConnection = MySQLdb.connect(host=mysqlhost,user=mysqluser,passwd=mysqlpasswd,db=stacksdb)
MyCursor = MyConnection.cursor()

lociinfo = {}

for locus in lociused:
	SQL = """SELECT readend,readpair,snps FROM catalog_index WHERE tag_id = %s"""%locus
	SQLLen = MyCursor.execute(SQL)
	OUT = MyCursor.fetchall()
	
	lociinfo[locus] = OUT
	
snps1 = []
snps2 = []
pairlist = []
lociread1 = []
lociread2 = []
pairlist1 = []
pairlist2 = []

#For each locus, create lists of reads 1 or 2 for #snps, list of locus names, and list of paired locus id.
#Pairlist contains a list of all tag_ids for the paired loci.

for locus in lociinfo:
	if lociinfo[locus][0][0] == 1:
		snps1.append(lociinfo[locus][0][2])
		pairlist.append(str(lociinfo[locus][0][1]))
		lociread1.append(locus)
		pairlist1.append(str(lociinfo[locus][0][1]))
	elif lociinfo[locus][0][0] == 2:
		snps2.append(lociinfo[locus][0][2])
		pairlist.append(str(lociinfo[locus][0][1]))
		lociread2.append(locus)
		pairlist2.append(str(lociinfo[locus][0][1]))
	else:
		pairlist.append(0)
		print "Locus %s not assinged to a read end"%locus
	

print "%d loci in the catalog"%len(lociused)
print "%d loci and %d snps from end 1"%(len(snps1),sum(snps1))
print "%d loci and %d snps from end 2"%(len(snps2),sum(snps2))

#Create a list of pairs in duple form, then go through and only keep those correctly assigned to the right read.
pairduples1 = list(zip(lociread1,pairlist1))
pairduples2 = list(zip(pairlist2,lociread2))


allpairduples = pairduples1
for pair in pairduples2:
	if pair in allpairduples:
		pass
	else:
		allpairduples.append(pair)
		

print "%d total fragments in catalog"%len(allpairduples)

#List of loci that have a read2 pair locus in the set of loci from read 2 in the catalog
matches1 = set(pairlist1) & set(lociread2)

print "%d fragments have snps on both read 1 and read 2"%len(matches1)

print "Sanity check...does this number,%d match how many loci in catalog?"%(len(allpairduples)+len(matches1))

outpairs = open("%s/MatchedReadPairs.csv"%outputdir,"w")

#Extract the pair druples for all fragments sequenced on both ends.
matchpairs = []
for locus in matches1:
	pos = [i for i, v in enumerate(allpairduples) if v[1]==locus]
	matchpairs.append(allpairduples[pos[0]])
	pairs = ",".join(allpairduples[pos[0]])
	outpairs.write("%s\n"%pairs)

vcf.close()
outpairs.close()

MyCursor.close()
MyConnection.close()


##########################################################################################
