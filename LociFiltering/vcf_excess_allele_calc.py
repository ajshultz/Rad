#! /usr/bin/env python

import re, sys, os, itertools, sets, getopt        # Load standard modules I often use
import MySQLdb   # Gives the ability to interact with the mysql database.

'''
This script will examine this distribution of loci in a stacks catalog that have individuals possessing more than 2 alleles, indicative of incorrect stacking.  In this case, only loci found in a given vcf file will be queried.
'''


def main(argv):
	try:
		opts,args = getopt.getopt(argv,'hv:l:u:p:o:s:m:',)
	except getopt.GetoptError:
		print "excess_allele_calc.py -v <Path to vcf file> -l <mysql host> -u <mysql user> -p <mysql password> -o <path to output directory (default .)> -s <Stacks database name> -m <minimum cut off depth for an allele to be considered for an individual (default = 0)>"
		sys.exit(2)
			
	mysqlhost = 'localhost'
	mysqluser = 'root'
	mysqlpasswd = ''
	stacksdb = ''
	outputdir = '.'
	depth = '0'
	vcffile = ''

	for opt, arg in opts:
		if opt == "-l":
			mysqlhost = arg
		elif opt == "-h":
			print "excess_allele_calc.py -v <Path to vcf file> -l <mysql host (default localhost)> -u <mysql user (default root)> -p <mysql password> -s <Stacks database name> -o <path to output directory (default .)> -m <minimum cut off depth for an allele to be considered for an individual (default = 0)>"
			sys.exit(2)
		elif opt == "-u":
			mysqluser = arg
		elif opt == "-p":
			mysqlpasswd = arg
		elif opt == "-o":
			outputdir = arg
		elif opt == "-s":
			stacksdb = arg
		elif opt == "-m":
			depth = arg
		elif opt == "-v":
			vcffile = arg


	morethan2 = open("%s/LociWithMoreThan2Alleles_mindepth%s_fromvcf.csv"%(outputdir,depth),"w")

	morethan2.write("catalog_id,#indivs>2,#indivs2,#indivs1\n")
	
	vcf = open(vcffile,"r")

	
	#Retrieve the catalog_ids that have more than 2 alleles in a single individual, incorporating the minimum stack depth cutoff
	MyConnection = MySQLdb.connect( host = mysqlhost, user = mysqluser, \
									passwd = mysqlpasswd, db = stacksdb)
	MyCursor = MyConnection.cursor()

	SQL = """SELECT DISTINCT catalog_id
FROM """+stacksdb+""".matches
WHERE catalog_id IN (SELECT *
     FROM (SELECT catalog_id
          FROM """+stacksdb+""".matches
          WHERE depth > """ +depth+"""
          GROUP BY sample_id,tag_id
          HAVING COUNT(tag_id) > 2)
     AS """+stacksdb+""");"""
	SQLLen = MyCursor.execute(SQL)  # returns the number of records retrieved

	mysql_output = MyCursor.fetchall()
	
	vcfloci = []

	for line in vcf:
		if line[0] != "#":
			newline = line.split("\t")
			if newline[2] not in vcfloci:
				vcfloci.append(newline[2])
			else:
				pass
		else:
			pass

	#Iterate through the catalog_ids retrieved above to retrive allele counts for all individuals that have that catalog_id.  Output counts of individuals with more than 2, 2, and 1 alleles to a csv file.

	for line in range(SQLLen):
		locus = str(mysql_output[line][0])
		if locus in vcfloci:
			locus_counts = [0,0,0]
			locus_sql = """SELECT COUNT(tag_id)
	FROM """+stacksdb+""".matches
	WHERE catalog_id = """+locus+"""
	GROUP BY sample_id,
		tag_id;"""
			locus_sql_len = MyCursor.execute(locus_sql)  # returns the number of records retrieve
			locus_sql_output = MyCursor.fetchall()
		
			for ind in range(locus_sql_len):
				if locus_sql_output[ind][0] == 1:
					locus_counts[2] += 1
				elif locus_sql_output[ind][0] == 2:
					locus_counts[1] += 1
				elif locus_sql_output[ind][0] > 2:
					locus_counts[0] += 1
			morethan2.write("%s,%d,%d,%d\n"%(locus,locus_counts[0],locus_counts[1],locus_counts[2]))

		else:
			pass
		
	
	morethan2.close()
	vcf.close()
	

if __name__ == "__main__":
	main(sys.argv[1:])