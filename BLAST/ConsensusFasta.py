#! /usr/bin/env python

import re, sys, os, itertools, sets, getopt        # Load standard modules I often use
import MySQLdb   # Gives the ability to interact with the mysql database.
import pickle
import vcf

"""
This script will take in a vcf file output from stacks, take the set of loci present in the vcf file, query the stacks database for the consensus sequence, and output a fasta file.
"""

#User input variables
vcffile = "/Users/allisonshultz/Dropbox/HFRad-Tags/HFPaired_Reduced/_PhyloPop_r.5p14/batch_1.vcf"
mysqlhost = "localhost"
mysqluser = "root"
mysqlpasswd = "hofi"
stacksdb = "HFdenovo_PairTestingReducedm3N3_radtags"
output = "../TestResults/testpy.fasta"


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

seqs = []

for locus in lociused:
	SQL = """SELECT seq FROM catalog_tags WHERE tag_id = %s"""%locus
	SQLLen = MyCursor.execute(SQL)
	OUT = MyCursor.fetchall()
	
	seqs.append(OUT[0][0])
	


#For each locus, create a list of all consensus sequences in the order of the locus list.


print "%d loci in the catalog"%len(lociused)


fastaout = open(output,"w")

#Write new fasta file.
for i in range(0,len(lociused)):	
	fastaout.write(">%s\n%s\n"%(lociused[i],seqs[i]))

vcf.close()
fastaout.close()

MyCursor.close()
MyConnection.close()
