#!/usr/bin/env python

import sys, getopt

#This is a script to check and see if there are any lines in the file of a different length than the input length. If found, those that are are written to an output file specified by the user. If none are found, a message stating this is printed to the screen and no output file is created.


def main(argv):
	try:
		opts,args = getopt.getopt(argv,'h1:l:o:',)
	except getopt.GetOptError:
		print "line_length_check.py -1 <Fastq file> -l <length> -o <output file name>"
		sys.exit(2)
			
	inputfile1 = ''
	length = 0
	outputfile = ''

	for opt, arg in opts:
		if opt == "-1":
			inputfile1 = arg
		elif opt == "-l":
			length = arg
		elif opt == "-o":
			outputfile = arg
		elif opt == "-h":
			print "line_length_check.py -1 <Fastq file> -l <length> -o <output file name>"
			sys.exit()


	end1 = open(inputfile1,"r")
	
	problem_loci = []

	linenum = 0
	
	for line in end1:
		linenum += 1
		if linenum % 4 == 1:
			name = line.strip('\n')
		elif linenum % 4 == 2:
			length1 = len(line.strip('\n'))
			if length1 != int(length):
				problem_loci.append([name,length1])
	
	if problem_loci == []:
		print "There are no sequences of a different length"
	else:
		out = open(outputfile,"w")
		for locus in problem_loci:
			out.write("%s\t%d\n"%(locus[0],locus[1]))
		out.close()
	
	end1.close()

if __name__ == "__main__":
		main(sys.argv[1:])