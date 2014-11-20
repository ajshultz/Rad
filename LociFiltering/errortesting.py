#! /usr/bin/env python

import re, sys, os, itertools, sets, getopt        # Load standard modules I often use
import MySQLdb   # Gives the ability to interact with the mysql database.
import indiv_read_pair_module
import pickle

mysqlhost = 'localhost'
mysqluser = 'root'
mysqlpasswd = 'hofi'
stacksdb = 'HF_FinalParamTesting_m4M1n1_radtags'
mincutoff = 4

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
	
print numsamples

for sample in range(1,numsamples+1):
	#Run the individual functions
	samres[sample] = indiv_read_pair_module.indreadassign(mysqlhost,mysqluser,mysqlpasswd,stacksdb,sample,mincutoff)
	samrespairs[sample] = indiv_read_pair_module.indlocuspairing(mysqlhost,mysqluser,mysqlpasswd,stacksdb,sample,mincutoff)

	#Write general numeric results to files
	print ('%d,%d,%d\n'%(sample,samres[sample][2],samres[sample][3]))


	#Find loci that have both read 1 and read 2 in an individual
	if samres[sample][1] == []:
		pass
	else:
		for loc in range(len(samres[sample][1])):
			samres[sample][1][loc] = str(samres[sample][1][loc])
		badloci = ','.join(samres[sample][1])
		print "Bad locus for sample %d, locus %s"%(sample,badloci)

	#Create a list of loci that incorrectly pair in an individual, and also write them to a file.
	if samrespairs[sample][2] == []:
		pass
	else:
		for locus in samrespairs[sample][2]:
			if locus not in mult_indivpairs:
				mult_indivpairs.append(locus)
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
				
	print catalog
	print "Sample %d is done!"%sample

