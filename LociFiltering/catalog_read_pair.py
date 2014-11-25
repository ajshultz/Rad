#! /usr/bin/env python

import re, sys, os, itertools, sets, getopt        # Load standard modules I often use
import MySQLdb   # Gives the ability to interact with the mysql database.
import indiv_read_pair_module
import pickle


def main(argv):
	try:
		opts,args = getopt.getopt(argv,'hl:u:p:o:m:s:',)
	except getopt.GetOptError:
		print "catalog_read_pair.py -l <mysql host> -u <mysql user> -p <mysql password> -o <path to output directory (default .)> -s <Stacks database name>  -m <minimum cut off depth for a stack to be included in catalogs>"
		sys.exit(2)
			
	mysqlhost = 'localhost'
	mysqluser = 'root'
	mysqlpasswd = ''
	stacksdb = ''
	mincutoff = 1
	outputdir = '.'

	for opt, arg in opts:
		if opt == "-l":
			mysqlhost = arg
		elif opt == "-h":
			print "catalog_read_pair.py -l <mysql host (default localhost)> -u <mysql user (default root)> -p <mysql password> -s <Stacks database name> -o <path to output directory (default .)> -m <minimum cut off depth for a stack to be included in catalogs (default 1)>"
			sys.exit(2)
		elif opt == "-u":
			mysqluser = arg
		elif opt == "-p":
			mysqlpasswd = arg
		elif opt == "-o":
			outputdir = arg
		elif opt == "-m":
			mincutoff = int(arg)
		elif opt == "-s":
			stacksdb = arg


	"""
	The first part of this script iterates through all individuals and runs them through indreadassign in order to assign all loci to reads, and map those read assignments back to the catalog tag numbers.  At the end, a file will be produced with the numbers of read 1 and read 2 snps for each sample, and any loci in each sample that are formed from both read 1 and read 2 loci (note that tag number is individual tag number).  In addition, a dictionary of catalog tag numbers with values that are lists of the reads for each individual is produced. Note that the catalog dictionary is saved as a pickle dump "catalogdict.p".
	"""

	samp_snps = open("%s/AllSampleSNPNumbers.csv"%outputdir,"w")
	overlap12 = open("%s/AllSampleReadOverlaps.csv"%outputdir,"w")

	samp_snps.write("Sample,Read1SNPS,Read2SNPS\n")
	overlap12.write("Sample,ProblemLoci\n")
	
	#Retrieve the number of samples from the mysql database
	MyConnection = MySQLdb.connect( host = mysqlhost, user = mysqluser, \
									passwd = mysqlpasswd, db = stacksdb)
	MyCursor = MyConnection.cursor()

	SQL = """SELECT file from samples;"""
	SQLLen = MyCursor.execute(SQL)  # returns the number of records retrieved
	
	numsamples = SQLLen
	

	samres = list(range(numsamples+1))
	samrespairs = list(range(numsamples+1))
	catalog = {}
	catalogpairs = {}

	mult_indivpairs = []

	#Files to write the locus pair results
	samp_pairinginfo = open("%s/AllSamplePairingInfo.csv"%outputdir,"w")
	samp_pairinginfo.write("Sample,NumberLoci,SinglyPaired,MultiplyPaired\n")

	mult_indivpairfile = open("%s/MultipleLociWithinInds.csv"%outputdir,"w")



	for sample in range(1,numsamples+1):
	#for sample in (20,21,47,48,53):
		#Run the individual functions
		samres[sample] = indiv_read_pair_module.indreadassign(mysqlhost,mysqluser,mysqlpasswd,stacksdb,sample,mincutoff)
		samrespairs[sample] = indiv_read_pair_module.indlocuspairing(mysqlhost,mysqluser,mysqlpasswd,stacksdb,sample,mincutoff)
	
		#Write general numeric results to files
		samp_snps.write('%d,%d,%d\n'%(sample,samres[sample][2],samres[sample][3]))
		pairres_write = (sample,len(samrespairs[sample][3]),len(samrespairs[sample][1]),len(samrespairs[sample][2]))
		samp_pairinginfo.write('%d,%d,%d,%d\n'%pairres_write)

		#Find loci that have both read 1 and read 2 in an individual
		if samres[sample][1] == []:
			pass
		else:
			for loc in range(len(samres[sample][1])):
				samres[sample][1][loc] = str(samres[sample][1][loc])
			badloci = ','.join(samres[sample][1])
			print "Bad locus for sample %d, locus %s"%(sample,badloci)
			overlap12.write("%d,%s,\n"%(sample,badloci))
	
		#Create a list of loci that incorrectly pair in an individual, and also write them to a file.
		if samrespairs[sample][2] == []:
			pass
		else:
			for locus in samrespairs[sample][2]:
				if locus not in mult_indivpairs:
					mult_indivpairs.append(locus)
					mult_indivpairfile.write('%s\n'%locus)
				else:
					pass
	
		#Create a catalog of read 1 or read 2 designations for catalog tag_ids for all individuals
		for locus in samres[sample][0]:
			if locus not in catalog:
				catalog[locus]=[samres[sample][0][locus]]
			else:
				catalog[locus].append(samres[sample][0][locus])

		#Create a catalog of all locus pairs for all individuals.  If a locus is found in more than one individual, the pairing is appended to the value for that key (the locus)
		for locus in samrespairs[sample][0]:
			if locus not in catalogpairs:
				catalogpairs[locus]=[samrespairs[sample][0][locus]]
			else:
				catalogpairs[locus].append(samrespairs[sample][0][locus])
	
		print "Sample %d is done!"%sample

	pickle.dump(catalog,open("%s/catalogdict.p"%outputdir,"w"))
	pickle.dump(catalogpairs,open("%s/catalogpairsdict.p"%outputdir,"w"))
	pickle.dump(mult_indivpairs,open("%s/list_multindivpairs.p"%outputdir,"w"))

	samp_snps.close()
	overlap12.close()
	samp_pairinginfo.close()
	mult_indivpairfile.close()

	##########################################################################################
	#Create a catalog of those loci with single read designations (collapsing them to a single value).

	singlecatalog = []

	for key in catalog:
		if sum(catalog[key])/float(len(catalog[key])) == 1:
			one = (1,key)
			singlecatalog.append(one)
		elif sum(catalog[key])/float(len(catalog[key])) == 2:
			two = (2,key)
			singlecatalog.append(two)
		else:
			three = (3,key)
			singlecatalog.append(three)
			print 'Locus %d has both read 1 and 2'%key

	pickle.dump(singlecatalog,open("%s/single_val_catalog.p"%outputdir,"w"))

	##########################################################################################

	"""With the paired read dictionary, this script will now iterate through to compress the dictionary into a single pair match, and identify those that have more than one locus match.
	"""

	singlecatalogpairs = {}
	catalogmultiplepairs = {}

	problemloci = open("%s/MultipleLociAcrossInds.csv"%outputdir,"w")
	problemloci.write("Catalogtag_id,,AllPairs\n")

	for key in catalogpairs:
		if sum(catalogpairs[key])%float(len(catalogpairs[key])) != 0:
			catalogmultiplepairs[key] = catalogpairs[key]
		else:
			singlecatalogpairs[key] = sum(catalogpairs[key])/len(catalogpairs[key])

	pickle.dump(singlecatalogpairs,open("%s/catalogsinglepairsdict.p"%outputdir,"w"))
	pickle.dump(catalogmultiplepairs,open("%s/catalogmultiplepairs.p"%outputdir,"w"))

	if catalogmultiplepairs != {}:
		print "Something is in catalogmultiplepairs"
		for key in catalogmultiplepairs:
			for loc in range(len(catalogmultiplepairs[key])):
				catalogmultiplepairs[key][loc] = str(catalogmultiplepairs[key][loc])
			loci = ",".join(catalogmultiplepairs[key])
			problemloci.write(("%d,,%s\n")%(key,loci))
	
	problemloci.close()


if __name__ == "__main__":
		main(sys.argv[1:])