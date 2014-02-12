#! /usr/bin/env python

import re, sys, os, itertools, sets, getopt        # Load standard modules I often use
import MySQLdb   # Gives the ability to interact with the mysql database.
import pickle
import vcf
import Pairwise_linkage_disequilibrium

"""
The first part of this script will obtain a list of locus catalog tag_id numbers from an inputted vcf file. It will match up loci to a stacks database, and pull information on which read and a paired locus (requires catalog_read_pair.py and mysql_database_pair_update.py having been run previously).  It will then output a number of stats, including how many loci are from read 1 or read 2, how many are matched to the same fragment, and a csv file of the matched loci where both ends are present in the library.
"""

#User input variables
vcffile = "/Users/allisonshultz/Dropbox/HFRad-Tags/HFPaired_Reduced/_PhyloPop_r.5p14/batch_1.vcf"
plinkped = "/Users/allisonshultz/Dropbox/HFRad-Tags/HFPaired_Reduced/_PhyloPop_r.5p14_allsitesoutput/batch_1.plink.ped"
plinkmap = "/Users/allisonshultz/Dropbox/HFRad-Tags/HFPaired_Reduced/_PhyloPop_r.5p14_allsitesoutput/batch_1.plink.map"
mysqlhost = "localhost"
mysqluser = "root"
mysqlpasswd = "hofi"
stacksdb = "HFdenovo_PairTestingReducedm3N3_radtags"
numsamples=78
outputdir = "TestResults"
popmap = "/Users/allisonshultz/Stacks/PopMaps_m/PopMapBasicPhylo_noMD.txt"



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
'''
This part of the script will calculate LD between pairs of SNPS.  LD will be calculated between all pairs of loci for each of the populations given by the popmap.
'''

#The PLINK ped file contains all genotype information for individuals in rows.  Note that the first 6 columns of each row are special.  Relevant here is the populations ID (first column) and individual ID (second column).  The PLINK map file contains locus information.  Most important is the second column, which has the locus catalog_id "_" position on the read.
ped = open(plinkped,"r")
map = open(plinkmap,"r")

#Get locus position information
locusid = []
catid = []
tagpos = []

for row in map:
	row = row.split("\t")
	locusid.append(row[1])
	loc = row[1]
	sep = loc.split("_")
	catid.append(sep[0])
	tagpos.append(sep[1])
	
#Get individual genotypes (make a dictionary), and make a dictionary of populations pop[list of individuals]
indgendict = {}
indpop = {}

for row in ped:
	row = row.strip()
	row = row.split("\t")
	if row[0] not in indpop:
		indpop[row[0]]=[row[1]]
	else:
		indpop[row[0]].append(row[1])
	indgendict[row[1]]=row[6:]
	numlets = len(row[6:])


#Make two lists, var1 with genotype 1 for a locus and var2 with genotype2 for a  locus.  The lists are in the same order as the catalog_ids and tag_postions from catid and tagpos.  The individual calls are in the order of indorder, which is extracted as the keys of the dictionary.  

index = range(0,numlets-1,2)

var1 = range(0,len(index))
var2 = range(0,len(index))

indorder = indgendict.keys()

for ind in indorder:
	for i in range(0,len(index)):
		if type(var1[i]) == list:
			var1[i].append(indgendict[ind][index[i]])
		else:
			var1[i] = [indgendict[ind][index[i]]]
		if type(var2[i]) == list:
			var2[i].append(indgendict[ind][index[i]+1])
		else:
			var2[i] = [indgendict[ind][index[i]+1]]		

genotypes = range(0,len(var1))

for i in range(0,len(var1)):
	genotypes[i] = tuple(zip(var1[i],var2[i]))
	
#Now, open the popmap and assign individuals into each of the populations listed in the popdict dictionary

pops = open(popmap,"r")


popdict = {}
for line in pops:
	line = line.strip()
	line = line.split("\t")
	pop = line[1]
	ind = line[0]
	if pop not in popdict:
		popdict[pop] = [ind]
	else:
		popdict[pop].append(ind)

#Extract the genotypes only for the individuals in each population, add them to the popgenotypes dictionary.

popgenotypes = {}

poppos = {}

for pop in popdict:
	popposlist = []
	for ind in popdict[pop]:
		indpos = [int(i) for i,v in enumerate(indorder) if v==ind]
		indpos = indpos[0]
		popposlist.append(indpos)
	poppos[pop]=popposlist

for pop in poppos:
	genset = []
	for locus in genotypes:
		newset = []
		for indpos in poppos[pop]:
			newset.append(locus[indpos])
		genset.append(tuple(newset))
	popgenotypes[pop]=genset


#Calculate LD for all pairs of SNPs for each population.  The resulting matrix is output into a separate file for each population (named by the population number) as a csv file. Note that this took about an hour for each population, so I will keep it commented out unless necessary. Note that loci are in same order as in PLINK files, so catalog and tag_id numbers can be obtained from Map file. 
# 
# popnames=popgenotypes.keys()
# 
# for popnum in range(0,len(popnames)):
# 	resfilename = (outputdir,"/pop_",popnames[popnum],".csv")
# 	resfilename = "".join(resfilename)
# 	resfile = open(resfilename,"w")
# 	for i in range(0,len(popgenotypes[popnames[popnum]])):
# 		liner2=[]
# 		for j in range(0,len(popgenotypes[popnames[popnum]])):
# 			try:
# 				r2 = Pairwise_linkage_disequilibrium.Pairwise_linkage_disequilibrium(popgenotypes[popnames[popnum]][i],popgenotypes[popnames[popnum]][j])
# 				r2 = round(r2['R_sq'],5)
# 			except:
# 				r2 = 9
# 			liner2.append(r2)
# 		liner2 = ','.join(str(k) for k in liner2)
# 		resfile.write("%s\n"%liner2)
# 	resfile.close()
	




ped.close()
map.close()
pops.close()


