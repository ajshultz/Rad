#! /usr/bin/env python

import re, sys, os, itertools, sets, getopt        # Load standard modules I often use
import MySQLdb   # Gives the ability to interact with the mysql database.
import pickle
import vcf

"""
This script will take in a vcf file output from stacks, take the set of loci present in the vcf file, query the stacks database for the consensus sequence, and output a fasta file.
"""
def main(argv):
	try:
		opts,args = getopt.getopt(argv,'hl:u:p:o:s:v:',)
	except getopt.GetOptError:
		print "ConsensusFasta.py -l <mysql host> -u <mysql user> -p <mysql password> -s <Stacks database name> -v <path to input vcf file> -o <path to output directory and filename (default ./vcf.fasta)> "
		sys.exit(2)
			
	mysqlhost = 'localhost'
	mysqluser = 'root'
	mysqlpasswd = ''
	stacksdb = ''
	output = './vcf.fasta'
	vcffile = ''

	for opt, arg in opts:
		if opt == "-l":
			mysqlhost = arg
		elif opt == "-h":
			print "ConsensusFasta.py -l <mysql host (default localhost)> -u <mysql user (default root)> -p <mysql password> -s <Stacks database name> -v <path to input vcf file> -o <path to output directory and filename (default ./vcf.fasta)> "
			sys.exit(2)
		elif opt == "-u":
			mysqluser = arg
		elif opt == "-p":
			mysqlpasswd = arg
		elif opt == "-o":
			output = arg
		elif opt == "-s":
			stacksdb = arg
		elif opt == "-v":
			vcffile = arg


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


if __name__ == "__main__":
		main(sys.argv[1:])