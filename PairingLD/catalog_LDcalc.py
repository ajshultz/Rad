#! /usr/bin/env python

import re, sys, os, itertools, sets, getopt        # Load standard modules I often use
import pickle
import Pairwise_linkage_disequilibrium


'''
This script will calculate LD between pairs of SNPS.  LD will be calculated between all pairs of loci for all individual in each of the populations given by the popmap, LD is calculated with the Pairwise_linkage_disequilibrium.py script from pypedia, obtained here: http://www.pypedia.com/index.php/Pairwise_linkage_disequilibrium.  Note r^2 is the LD metric output by this script.  Required input is an output directory to write results, popmap file, plink ped file, and plink map file.
'''

def main(argv):
	try:
		opts,args = getopt.getopt(argv,'ho:m:P:M:',)
	except getopt.GetOptError:
		print "catalog_read_pair.py -o <path to output directory (default .)> -m <Popmap file> -P <Plink ped file> -M <Plink map file>"
		sys.exit(2)
			
	outputdir = '.'
	plinkped = ''
	plinkmap = ''
	popmap = ''

	for opt, arg in opts:
		if opt == "-h":
			print "catalog_read_pair.py -o <path to output directory (default .)> -m <Popmap file> -P <Plink ped file> -M <Plink map file>"
			sys.exit(2)
		elif opt == "-o":
			outputdir = arg
		elif opt == "-m":
			popmap = arg
		elif opt == "-P":
			plinkped = arg
		elif opt == "-M":
			plinkmap = arg



	#The PLINK ped file contains all genotype information for individuals in rows.  Note that the first 6 columns of each row are special.  Relevant here is the populations ID (first column) and individual ID (second column).  The PLINK map file contains locus information.  Most important is the second column, which has the locus catalog_id "_" position on the read.
	ped = open(plinkped,"r")
	map = open(plinkmap,"r")

	#Get locus position information
	locusid = []
	catid = []
	tagpos = []

	for row in map:
		if row[0] == "#":
			pass
		else:
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
		if row[0] == "#":
			pass
		else:
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


	#Calculate LD for all pairs of SNPs for each population.  The resulting matrix is output into a separate file for each population (named by the population number) as a csv file. Note that this took about an hour for each population, so I will keep it commented out unless necessary. Note that loci are in same order as in PLINK files, so catalog and tag_id numbers can be obtained from Map file. Note that there are some cases where r2 cannot be calculated due to a division by 0.  In these cases, a 9 is given for the r2 value as a place holder. 

	popnames=popgenotypes.keys()

	for popnum in range(0,len(popnames)):
		resfilename = (outputdir,"/pop_",popnames[popnum],".csv")
		resfilename = "".join(resfilename)
		resfile = open(resfilename,"w")
		for i in range(0,len(popgenotypes[popnames[popnum]])):
			liner2=[]
			for j in range(0,len(popgenotypes[popnames[popnum]])):
				try:
					r2 = Pairwise_linkage_disequilibrium.Pairwise_linkage_disequilibrium(popgenotypes[popnames[popnum]][i],popgenotypes[popnames[popnum]][j])
					r2 = round(r2['R_sq'],5)
				except:
					r2 = 9
				liner2.append(r2)
			liner2 = ','.join(str(k) for k in liner2)
			resfile.write("%s\n"%liner2)
		resfile.close()
	




	ped.close()
	map.close()
	pops.close()


if __name__ == "__main__":
		main(sys.argv[1:])