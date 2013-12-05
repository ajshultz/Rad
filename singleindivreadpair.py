#! /usr/bin/env python

""" 
indivreadpair.py
using the mysql stacks database, with its table 'unique_tags',
"""

import re        # Load regular expression module
import MySQLdb   # must be installed separately
import sets
import sys

stacksdb = "HFdenovo_PairTesting_radtags"
samp_id = "20"	

outfile="~/Stacks/HF_PairedTestig_indiv1.txt"


# Create the database connection
# Often you will want to use a variable instead of a fixed string 
# for the database name
MyConnection = MySQLdb.connect( host = "localhost", user = "root", \
                                passwd = "hofi", db = stacksdb)
MyCursor = MyConnection.cursor()

SQL = """SELECT seq_id,tag_id from unique_tags WHERE sample_id=%s;"""%samp_id
SQLLen = MyCursor.execute(SQL)  # returns the number of records retrieved
#print SQLLen

# MyCursor is now "loaded" with the results of the SQL command
# AllOut will become a list of all the records selected
AllOut = MyCursor.fetchall()   
#print AllOut

read1 = []
read1keys = []
read1values = []
read2 = []
read2keys=[]
read2values=[]

#read1loci and read2loci will contain each unique locus (not repeated).
read1loci = []
read2loci = []

#Go through all reads from an individual and assign them based on whether they are from read1 or read2 to a list of the sequence id and tag id for that individual, a list of just the sequence id and a list of just the tag id..these are then zipped into a dictionary with the key as the seq_id and the value as the tag_id

for index in range(SQLLen):
	if AllOut[index][0] == "":
		pass
	elif AllOut[index][0][-1] == "1":
		read1.append(AllOut[index])
		read1keys.append(AllOut[index][0][:-2])
		read1values.append(AllOut[index][1])
		if AllOut[index][1] in read1loci:
			pass
		else:
			read1loci.append(AllOut[index][1])
	elif AllOut[index][0][-1] == "2":
		read2.append(AllOut[index])
		read2keys.append(AllOut[index][0][:-2])
		read2values.append(AllOut[index][1])
		if AllOut[index][1] in read2loci:
			pass
		else:
			read2loci.append(AllOut[index][1])
		
	else:
		pass


read1dict = dict(zip(read1keys,list(read1values)))
read2dict = dict(zip(read2keys,read2values))

#Identify any tags that were formed from reads 1 and read 2.

#print "Loci formed from read1 and read2: ",list(set(read1values) & set(read2values)),"\n\n\n"


#print "Read1: ",read1loci,"\n\n\n"
#print "Read2: ",read2loci,"\n\n\n"


#Now we are going to pull data from the tag_index table, which contains the individual tag_id, catalog_id (catalog tag_id), sequencing depth, and snps. We only pull the records for a particular individual.

SQL2 = """SELECT tag_id,catalog_id,depth,snps from tag_index WHERE sample_id=%s;"""%samp_id
SQLLen2 = MyCursor.execute(SQL2)

AllOut2 = MyCursor.fetchall()

#First, with those records we pull how many snps were in read 1 and read 2. At the same time, we create a dictionary of catalog_id with which read it was from.

snps1 = int(0)
snps2 = int(0)
#cat_read = []
indloci=[]

for index in range(SQLLen2):
	indloci.append([AllOut2[index][0]])
	if AllOut2[index][0] in read1loci:
		snps1 = snps1+AllOut2[index][3]
#		cat_read[AllOut2[index][1]]=1
	elif AllOut2[index][0] in read2loci:
		snps2 = snps2+AllOut2[index][3]
#		cat_read[AllOut2[index][1]]=2
	else:
		print "Something is wrong with locus %d.\n"%AllOut2[index][0]
		
#print "Read 1 has %d Snps.\n\n\n"%snps1
#print "Read 2 has %d Snps.\n\n\n"%snps2	


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
	

allkeys = allreads.keys()

for key in allkeys:
	if len(allreads[key]) < 2:
		del allreads[key]
	else:
		pass

#print allreads
unique = {}

for key in allreads:
	if allreads[key] not in unique:
		unique	[allreads[key]]=1
	else:
		unique[allreads[key]]+=1


uniquepairs = unique.keys()

#print [i for i, v in enumerate(uniquepairs) if v[0]==1 or v[1]==1]

singlepairs = []
multiplepairs = []
cat_pairindex = {}

for loc in range(len(indloci)):
	pos1 = [i for i, v in enumerate(uniquepairs) if v[0]==indloci[loc][0]]
	pos2 = [i for i, v in enumerate(uniquepairs) if v[1]==indloci[loc][0]]
	if len(pos1)+len(pos2) == 1:
		singlepairs.append(AllOut2[indloci[loc][0]-1[1])
		if len(pos1) == 1:
			cat_pairindex[AllOut2[indloci[loc][0]-1][1]] = AllOut2[uniquepairs[pos1[0]][1]-1][1]
		elif len(pos2) == 1:
			cat_pairindex[AllOut2[indloci[loc][0]-1][1]] = AllOut2[uniquepairs[pos2[0]][0]-1][1]
		else:
			print "Something is wrong"
	elif len(pos1)+len(pos2) > 1:
		multiplepairs.append(AllOut2[indloci[loc][0]-1][1])
	
print AllOut2

print cat_pairindex
print multiplepairs
print singlepairs
print "There are %d loci with only one partner.\n"%len(singlepairs)
print "There are %d loci with more than one partner.\n"%len(multiplepairs)




MyCursor.close()
MyConnection.close()