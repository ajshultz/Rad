#!/usr/bin/env python

import sys, getopt

#This is a script to trim all sequences in a fastq file to a specified length.


def main(argv):
	try:
		opts,args = getopt.getopt(argv,'h1:l:o:',)
	except getopt.GetOptError:
		print "line_length_trimming.py -1 <input Fastq file> -l <length> -o <output Fastq file name>"
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
			print "line_length_trimming.py -1 <input Fastq file> -l <length> -o <output Fastq file name>"
			sys.exit()
		

	fastq = open(inputfile1,"r")
	out = open(outputfile,"w")

	linenum = 0
	
	for line in fastq:
		linenum += 1
		if linenum % 4 == 1:
			out.write(line)
		elif linenum % 4 == 2:
			newline = line.strip('\n')
			length1 = len(newline)
			if length1 == int(length):
				out.write("%s\n"%(newline))
			else:
				diff = length1-int(length)
				if diff <= 0:
					print "There is a sequence shorter than the length given, it is only %d bp long"%length1
				else:
					newline = newline[:-diff]
					out.write("%s\n"%(newline))
		elif linenum % 4 == 3:
			out.write(line)
		elif linenum % 4 == 0:
			newline = line.strip('\n')
			length1 = len(newline)
			if length1 == int(length):
				out.write("%s\n"%(newline))
			else:
				diff = length1-int(length)
				if diff <= 0:
					print "There is a sequence shorter than the length given, it is only %d bp long"%length1
				else:
					newline = newline[:-diff]
					out.write("%s\n"%(newline))
			

	out.close()
	fastq.close()

if __name__ == "__main__":
		main(sys.argv[1:])