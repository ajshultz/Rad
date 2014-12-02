#! /usr/bin/env python

import re, sys, os, itertools, sets, getopt        # Load standard modules I often use

'''
This script will concatenate the results of multiple parameter analyses and produce csv files that can be read into other statistical programs for graphical display.  Note that the same samples should have been used for each stacks run.
'''


def main(argv):
	try:
		opts,args = getopt.getopt(argv,'hi:o:m:',)
	except getopt.GetoptError:
		print "all_param_results.py -i <input dir containing different parameter result directories> -o <path to output directory (default .)> -m <minimum cut off depth used in previous scripts (default = 0)>"
		sys.exit(2)
	inputdir = '.'
	outputdir = '.'
	depth = '0'

	for opt, arg in opts:
		if opt == "-h":
			print "all_param_results.py -i <input dir containing different parameter result directories> -o <path to output directory (default .)> -m <minimum cut off depth used in previous scripts (default = 0)>"
			sys.exit(2)
		elif opt == "-i":
			inputdir = arg
		elif opt == "-o":
			outputdir = arg
		elif opt == "-m":
			depth = arg

	subdirs = next(os.walk(inputdir))[1]

	#Create a dictionary to hold a list of results from all subdirectories for each parameter, and an empty list for each parameter of interes
	param_res = {}
	
	param_res["sample"] = []
	param_res["num_loci"] = []
	param_res["num_single_pair_loci"] = []
	param_res["num_mult_pair_loci"] = []
	param_res["num_excess_alleles_loci"] = []

	for dir in subdirs:
		subdir_path = inputdir+"/"+dir
		pairinginfo = open(subdir_path+"/AllSamplePairingInfo.csv","r")
		excessinfo = open(subdir_path+"/LociWithMoreThan2Alleles_PerSample_mindepth"+depth+".csv","r")
		
		#Create empty lists to hold the results of each particular file
		samples = []
		num_loci = []
		num_singles = []
		num_mult = []
		num_excess = []
		
		#Add all individual results to lists
		for line in pairinginfo:
			line = line.strip().split(",")
			if line[0] == "Sample":
				pass
			else:
				samples.append(line[0])
				num_loci.append(line[1])
				num_singles.append(line[2])
				num_mult.append(line[3])
		
		for line in excessinfo:
			line = line.strip().split(",")
			if line[0] == "Sample":
				pass
			else:
				num_excess.append(line[1])
		
		#Add result lists back to dictionary
		param_res["sample"].append(samples)
		param_res["num_loci"].append(num_loci)
		param_res["num_single_pair_loci"].append(num_singles)
		param_res["num_mult_pair_loci"].append(num_mult)
		param_res["num_excess_alleles_loci"].append(num_excess)		
		
		pairinginfo.close()
		excessinfo.close()

	#Quick check to make sure each file has the same number of samples.
	num_samples = []
	for i in range(len(subdirs)):
		num_samples.append(len(param_res["sample"][i]))
	if sum(num_samples)%len(num_samples) != 0:
		sys.exit("Oh, no! The number of samples doesn't match among files!")	
	else:
		pass
	
	#Create a single list of sample_ids
	samples = param_res["sample"][0]


	num_loci_out = open(outputdir+"/Num_Loci.csv","w")
	num_singles_out = open(outputdir+"/Num_Single_Pair_Loci.csv","w")
	num_mult_out = open(outputdir+"/Num_Multiple_Pair_Loci.csv","w")
	num_excess_out = open(outputdir+"/Num_Excess_Allele_Loci.csv","w")
	per_singles_out = open(outputdir+"/Percentage_Single_Pair_Loci.csv","w")
	per_mult_out = open(outputdir+"/Percentage_Multiple_Pair_Loci.csv","w")
	per_excess_out = open(outputdir+"/Percentage_Excess_Allele_Loci.csv","w")

	subdir_header = "sample,"+",".join(subdirs)+"\n"
	
	num_loci_out.write(subdir_header)
	num_singles_out.write(subdir_header)
	num_mult_out.write(subdir_header)
	num_excess_out.write(subdir_header)
	per_singles_out.write(subdir_header)
	per_mult_out.write(subdir_header)
	per_excess_out.write(subdir_header)

	
	ind_num_loci = []
	ind_num_singles = []
	ind_num_mult = []
	ind_num_excess = []
	ind_per_singles = []
	ind_per_mult = []
	ind_per_excess = []
	
	for ind in range(len(samples)):
		ind_num_loci.append([])
		ind_num_singles.append([])
		ind_num_mult.append([])
		ind_num_excess.append([])
		ind_per_singles.append([])
		ind_per_mult.append([])
		ind_per_excess.append([])
		
		for i in range(len(subdirs)):
			ind_num_loci[ind].append(param_res["num_loci"][i][ind])
			ind_num_singles[ind].append(param_res["num_single_pair_loci"][i][ind])
			ind_num_mult[ind].append(param_res["num_mult_pair_loci"][i][ind])
			ind_num_excess[ind].append(param_res["num_excess_alleles_loci"][i][ind])
			ind_per_singles[ind].append(str((float(param_res["num_single_pair_loci"][i][ind])/float(param_res["num_loci"][i][ind]))))
			ind_per_mult[ind].append(str((float(param_res["num_mult_pair_loci"][i][ind])/float(param_res["num_loci"][i][ind]))))
			ind_per_excess[ind].append(str((float(param_res["num_excess_alleles_loci"][i][ind])/float(param_res["num_loci"][i][ind]))))
	
	for ind in range(len(samples)):
		num_loci_out.write("%s,%s\n"%(samples[ind],",".join(ind_num_loci[ind])))
		num_singles_out.write("%s,%s\n"%(samples[ind],",".join(ind_num_singles[ind])))
		num_mult_out.write("%s,%s\n"%(samples[ind],",".join(ind_num_mult[ind])))
		num_excess_out.write("%s,%s\n"%(samples[ind],",".join(ind_num_excess[ind])))
		per_singles_out.write("%s,%s\n"%(samples[ind],",".join(ind_per_singles[ind])))
		per_mult_out.write("%s,%s\n"%(samples[ind],",".join(ind_per_mult[ind])))
		per_excess_out.write("%s,%s\n"%(samples[ind],",".join(ind_per_excess[ind])))

	
	num_loci_out.close()
	num_singles_out.close()
	num_mult_out.close()
	num_excess_out.close()
	per_singles_out.close()
	per_mult_out.close()
	per_excess_out.close()

	
	
	
if __name__ == "__main__":
	main(sys.argv[1:])