#!/usr/bin/env python

import sys, getopt

#This is a function to check and see if there are any lines in the file of a different length than the input length. Those that are are written to an output file specified by the user.


def main(argv):
	try:
		opts,args = getopt.getopt(argv,'h1:l:o:',)
	except getopt.GetOptError:
		print "LineLengthcheck.py -1 <Fastq file> -l <length> -o <output>"
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
		


	end1 = open(inputfile1,"r")
	out = open(outputfile,"w")

	linenum = 0
	
	for line in end1:
		linenum += 1
		if linenum % 4 == 1:
			name = line.strip('\n')
		elif linenum % 4 == 2:
			length1 = len(line.strip('\n'))
			if length1 != int(length):
				out.write("%s\t%d\n"%(name,length1))

		

	out.close()
	end1.close()

if __name__ == "__main__":
		main(sys.argv[1:])