#! /usr/bin/env python

""" 
indivreadpair.py
using the mysql stacks database, with the table 'unique_tags' and 'index_tags',
"""

import re, sys, os, itertools, sets        # Load standard modules I often use
import MySQLdb   # Gives the ability to interact with the mysql database.

"""
This function will access the mysql database for an individual given by samp_id, extract the read information from the unique_tags read id's as to whether the reads are from the first or second read, then get the information from index_tags to extract the number of snps for that locus, and make a dictionary with the catalog_id as the key and the read as the value.  The mincutoff value represents the minimum depth of a locus to be added to the catalog dictionary.  The function will return a list of the catalog_id[read] dictionary, list of overlapping loci (loci found in both read 1 and read 2), the number of snps in read 1 and the number of snps in read 2.
"""
def indreadassign(mysqlhost,mysqluser,mysqlpasswd,stacksdb,samp_id,mincutoff):

	# Create the database connection
	# Often you will want to use a variable instead of a fixed string 
	# for the database name
	MyConnection = MySQLdb.connect( host = mysqlhost, user = mysqluser, \
									passwd = mysqlpasswd, db = stacksdb)
	MyCursor = MyConnection.cursor()

	SQL = """SELECT seq_id,tag_id from unique_tags WHERE sample_id=%s;"""%samp_id
	SQLLen = MyCursor.execute(SQL)  # returns the number of records retrieved

	# MyCursor is now "loaded" with the results of the SQL command
	# AllOut will become a list of all the records selected
	AllOut = MyCursor.fetchall()   
	#print AllOut

	#read1loci and read2loci will contain each unique locus (not repeated).

	read1 = []
	read2 = []
	read1loci = []
	read2loci = []

	#Go through all reads from an individual and assign them based on whether they are from read1 or read2 to a list of the sequence id and tag id for that individual, a list of just the sequence id and a list of just the tag id..these are then zipped into a dictionary with the key as the seq_id and the value as the tag_id

	for index in range(SQLLen):
		if AllOut[index][0] == "":
			pass
		elif AllOut[index][0][-1] == "1":
			read1.append(AllOut[index])
			if AllOut[index][1] in read1loci:
				pass
			else:
				read1loci.append(AllOut[index][1])
		elif AllOut[index][0][-1] == "2":
			read2.append(AllOut[index])
			if AllOut[index][1] in read2loci:
				pass
			else:
				read2loci.append(AllOut[index][1])
		
		else:
			pass


	#Identify any tags that were formed from reads 1 and read 2.
	
	overlaploci = list(set(read1loci) & set(read2loci))

	#print "Loci formed from read1 and read2: ",overlaploci,"\n\n\n"
	#print "Read1: ",read1loci,"\n\n\n"
	#print "Read2: ",read2loci,"\n\n\n"


	#Now we are going to pull data from the tag_index table, which contains the individual tag_id, catalog_id (catalog tag_id), sequencing depth, and snps. We only pull the records for a particular individual.

	SQL2 = """SELECT tag_id,catalog_id,depth,snps from tag_index WHERE sample_id=%s;"""%samp_id
	SQLLen2 = MyCursor.execute(SQL2)

	AllOut2 = MyCursor.fetchall()

	#First, with those records we pull how many snps were in read 1 and read 2. At the same time, we create a dictionary of catalog_id with which read it was from.

	snps1 = int(0)
	snps2 = int(0)
	cat_read = {}

	for index in range(SQLLen2):
		if AllOut2[index][2] >= mincutoff:
			if AllOut2[index][0] in read1loci:
				snps1 = snps1+AllOut2[index][3]
				cat_read[AllOut2[index][1]]=1
			elif AllOut2[index][0] in read2loci:
				snps2 = snps2+AllOut2[index][3]
				cat_read[AllOut2[index][1]]=2
			else:
				print "Something is wrong with locus %d.\n"%AllOut2[index][0]
		else:
			pass
		
	#print "Read 1 has %d Snps.\n\n\n"%snps1
	#print "Read 2 has %d Snps.\n\n\n"%snps2	

	MyCursor.close()
	MyConnection.close()

	return [cat_read,overlaploci,snps1,snps2]




##########################################################################################
"""
This function will access the mysql database for an individual given by samp_id, extract the read information from the unique_tags read id's as to whether the reads are from the first or second read, then use the sequencer position name to match paired reads.  This function will then identify which loci have more than one match, indicating that they likely represent paralogous regions that have incorrectly grouped together.  The function will return a catalog of catalog tag_ids with their pair (if there was only one locus paired), a list of singly paired catalog tag_ids, a list of multiply paired catalog tag_ids, and a list of possible loci for that individual.
"""

def indlocuspairing(mysqlhost,mysqluser,mysqlpasswd,stacksdb,samp_id,mincutoff):
	# Create the database connection
	# Often you will want to use a variable instead of a fixed string 
	# for the database name
	MyConnection = MySQLdb.connect( host = mysqlhost, user = mysqluser, \
									passwd = mysqlpasswd, db = stacksdb)
	MyCursor = MyConnection.cursor()

	SQL = """SELECT seq_id,tag_id from unique_tags WHERE sample_id=%s;"""%samp_id
	SQLLen = MyCursor.execute(SQL)  # returns the number of records retrieved

	# MyCursor is now "loaded" with the results of the SQL command
	# AllOut will become a list of all the records selected
	AllOut = MyCursor.fetchall()   
	
	#Create a dictionary with Illumina read IDs as keys and individual tag_id as values.  This is done separately for read 1 and read 2.
	
	read1keys = []
	read1values = []
	read2keys=[]
	read2values=[]
	
	for index in range(SQLLen):
		if AllOut[index][0] == "":
			pass
		elif AllOut[index][0][-1] == "1":
			read1keys.append(AllOut[index][0][:-2])
			read1values.append(AllOut[index][1])

		elif AllOut[index][0][-1] == "2":
			read2keys.append(AllOut[index][0][:-2])
			read2values.append(AllOut[index][1])
		else:
			pass
	
	read1dict = dict(zip(read1keys,list(read1values)))
	read2dict = dict(zip(read2keys,read2values))

	#Pull data from the tag_index table, which contains the individual tag_id, catalog_id (catalog tag_id), sequencing depth, and snps. We only pull the records for a particular individual.

	SQL2 = """SELECT tag_id,catalog_id,depth,snps from tag_index WHERE sample_id=%s;"""%samp_id
	SQLLen2 = MyCursor.execute(SQL2)

	AllOut2 = MyCursor.fetchall()
	
	
	#Get list of all loci for the individual and create a dictionary of Mysql output for loci that pass the minimum sampling depth filter, with the individual tag_id as the key.
	indloci=[]
	mysql_tagindex={}
	for index in range(SQLLen2):
		if AllOut2[index][2] >= mincutoff:
			indloci.append([AllOut2[index][0]])
			mysql_tagindex[AllOut2[index][0]]=AllOut2[index]
		else:
			pass
			#print "Locus %d doesn't not meet minumum depth requirments, it is %d"%(index,AllOut2[index][2])


	#Make a dictionary of all possible reads (1 and 2), then turn each pair into a tuple.
	allreads={}
	for key in read1dict:
		allreads.setdefault(key,[]).append(read1dict[key])
	for key in read2dict:
		if key in allreads:
			allreads[key].append(read2dict[key])
		else:
			pass
	for key in allreads:
		allreads[key]=tuple(allreads[key])
	
	#Keep only those loci that are actually paired
	allkeys = allreads.keys()
	for key in allkeys:
		if len(allreads[key]) < 2:
			del allreads[key]
		else:
			pass

	#Create a dictionary of all unique pair combinations, the values for the dictionary are the number of reads for that pairing.
	unique = {}
	for key in allreads:
		if allreads[key] not in unique:
			unique	[allreads[key]]=1
		else:
			unique[allreads[key]]+=1

	#Get a list of all unique pairs
	uniquepairs = unique.keys()

	#Create a list of loci that have a single match, a list of loci that have more than one match, and a catalog of the loci with only one pair.  Note that the tag_id has been translated to the catalog_id at this point for later comparison among individuals.
	
	singlepairs = []
	multiplepairs = []
	cat_pairindex = {}

	#position 1 and position 2 have to be searched independently.
	for loc in range(len(indloci)):
		pos1 = [i for i, v in enumerate(uniquepairs) if v[0]==indloci[loc][0]]
		pos2 = [i for i, v in enumerate(uniquepairs) if v[1]==indloci[loc][0]]
		
		if len(pos1)+len(pos2) == 1:
			singlepairs.append(mysql_tagindex[indloci[loc][0]][1])
			if len(pos1) == 1:
				cat_pairindex[mysql_tagindex[indloci[loc][0]][1]] = mysql_tagindex[uniquepairs[pos1[0]][1]][1]
			elif len(pos2) == 1:
				cat_pairindex[mysql_tagindex[indloci[loc][0]][1]] = mysql_tagindex[uniquepairs[pos2[0]][0]][1]
			else:
				print "Something is wrong"
		elif len(pos1)+len(pos2) > 1:
			multiplepairs.append(mysql_tagindex[indloci[loc][0]][1])	
	return(cat_pairindex,singlepairs,multiplepairs,indloci)