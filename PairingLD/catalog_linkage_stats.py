#! /usr/bin/env python

import re, sys, os, itertools, sets, getopt        # Load standard modules I often use
import MySQLdb   # Gives the ability to interact with the mysql database.
import pickle
import vcf

"""
The first part of this script will obtain a list of locus catalog tag_id numbers from an inputted vcf file. It will match up loci to a stacks database, and pull information on which read and a paired locus (requires catalog_read_pair.py and mysql_database_pair_update.py having been run previously).  It will then output a number of stats in the file catalog_linkage_stats.csv, including how many loci are from read 1 or read 2, how many are matched to the same fragment, and a csv file of the matched loci where both ends are present in the library.
"""

def main(argv):
	try:
		opts,args = getopt.getopt(argv,'hl:u:s:p:o:v:',)
	except getopt.GetOptError:
		print "catalog_read_pair.py -l <mysql host> -u <mysql user> -p <mysql password> -s <Stacks database name> -o <path to output directory (default .)> -v <VCF file>"
		sys.exit(2)
			
	mysqlhost = 'localhost'
	mysqluser = 'root'
	mysqlpasswd = ''
	stacksdb = ''
	outputdir = '.'
	vcffile = ''

	for opt, arg in opts:
		if opt == "-l":
			mysqlhost = arg
		elif opt == "-h":
			print "catalog_read_pair.py -l <mysql host> -u <mysql user> -p <mysql password> -s <Stacks database name> -o <path to output directory (default .)> -v <VCF file>"
			sys.exit(2)
		elif opt == "-u":
			mysqluser = arg
		elif opt == "-p":
			mysqlpasswd = arg
		elif opt == "-s":
			stacksdb = arg
		elif opt == "-o":
			outputdir = arg
		elif opt == "-v":
			vcffile = arg

	vcf = open(vcffile,"r")
	outstats = open('%s/catalog_linkage_stats.csv'%(outputdir),"w")

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


	
	#Retrieve the number of samples from the mysql database
	MyConnection = MySQLdb.connect( host = mysqlhost, user = mysqluser, \
									passwd = mysqlpasswd, db = stacksdb)
	MyCursor = MyConnection.cursor()

	SQL = """SELECT file from samples;"""
	SQLLen = MyCursor.execute(SQL)  # returns the number of records retrieved

	numsamples = SQLLen



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
	

	outstats.write("total_catalog_loci,%d\n"%len(lociused))
	outstats.write("end_1_loci,%d\n"%(len(snps1)))
	outstats.write("end_2_lcoi,%d\n"%(len(snps2)))
	outstats.write("end_1_snps,%d\n"%(sum(snps1)))
	outstats.write("end_2_snps,%d\n"%(sum(snps2)))

	#Create a list of pairs in duple form, then go through and only keep those correctly assigned to the right read.
	pairduples1 = list(zip(lociread1,pairlist1))
	pairduples2 = list(zip(pairlist2,lociread2))


	allpairduples = pairduples1
	for pair in pairduples2:
		if pair in allpairduples:
			pass
		else:
			allpairduples.append(pair)
		

	outstats.write("total_catalog_fragments,%d\n"%len(allpairduples))

	#List of loci that have a read2 pair locus in the set of loci from read 2 in the catalog
	matches1 = set(pairlist1) & set(lociread2)

	outstats.write("fragments_with_read1_and_read2_snps,%d\n"%len(matches1))

	#print "Sanity check...does this number,%d match how many loci in catalog?"%(len(allpairduples)+len(matches1))

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
	outstats.close()

	MyCursor.close()
	MyConnection.close()


if __name__ == "__main__":
		main(sys.argv[1:])