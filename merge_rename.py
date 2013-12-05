#!/usr/bin/env python

import sys, getopt

#Note that this function is for taking the output of process_radtags from Stacks, and formatting paired-end reads from double-digest RADseq to be used for de novo assembly.  Input is read 1 and read 2 for an individual, and possible an extra file with unpaired read 1 data. This script will concatenate these files, including an option to re-add the read information to the illumina position info (_1 or _2), and to trim n bases off of read 2 to make it the same length as read 1.  


def main(argv):
	try:
		opts,args = getopt.getopt(argv,'h1:2:s:o:n:t:',)
	except getopt.GetOptError:
		print "merge_rename.py -1 <End 1 Fastq Input File> -2 <End 2 Fastq Input File> -s <Unpaired Fastq Input File> -o <Output File> -n <T or F, Add _1 or _2 back to read names> -t < T or F, trim n given bases off end of second read>"
		sys.exit(2)
			
	inputfile1 = ''
	inputfile2 = ''
	inputfilesingles = ''
	outputfile = ''
	rename = 'F'
	trim = 'F'

	for opt, arg in opts:
		if opt == "-1":
			inputfile1 = arg
		elif opt == "h":
			print "merge_rename.py -1 <End 1 Input File> -2 <End 2 Input File> -s <Unpaired Input Files> -o <Output File> -n <Add :1 read names> -t <Trim given bases off end of second read>"
		elif opt == "-2":
			inputfile2 = arg
		elif opt == "-s":
			inputfilesingles = arg
		elif opt == "-o":
			outputfile = arg
		elif opt == "-n":
			rename = arg
		elif opt == "-t":
			trim = arg
		
	print "Input file 1 is ",inputfile1

	end1 = open(inputfile1,"r")
	end2 = open(inputfile2,"r")
	if inputfilesingles != '':
		singles = open(inputfilesingles,"r")
	else:
		pass
	out = open(outputfile,"w")

	linenum = 0
	
	for line in end1:
		linenum += 1
		if linenum % 4 == 1:
			if rename == "T":
				newline = line.strip('\n')
				elements = newline.split('_')
				elements[-1] = "1"
				newline = '_'.join(elements)
				out.write(newline+'\n')
			else:
				out.write(line)
		elif linenum % 4 == 2:
			length1 = len(line.strip('\n'))
			out.write(line)
		elif linenum % 4 == 3:
			out.write(line)
		elif linenum % 4 == 0:
			out.write(line)

	linenum = 0

	for line in end2:
		linenum += 1
		if linenum % 4 == 1:
			if rename == "T":
				newline = line.strip('\n')
				elements = newline.split('_')
				elements[-1] = "2"
				newline = '_'.join(elements)
				out.write(newline+'\n')
			else:
				out.write(line)
		if linenum % 4 == 2:
			if trim == "T":
				newline = line.strip('\n')
				length2 = len(newline)
				diff = length2 - length1
				newline = newline[:-diff]
				out.write(newline + '\n')
			else:
				out.write(line)
		if linenum % 4 == 3:
			out.write(line)
		if linenum % 4 == 0:
			if trim == "T":
				newline = line.strip('\n')
				length2 = len(newline)
				diff = length2 - length1
				newline = newline[:-diff]
				out.write(newline + '\n')
			else:
				out.write(line)
	
	print inputfilesingles
	if inputfilesingles != '':
	
		linenum = 0
	
		for line in singles:
			linenum += 1
			if linenum % 4 == 1:
				if rename == "T":
					newline = line.strip('\n')
					elements = newline.split('_')
					elements[-1] = "1"
					newline = '_'.join(elements)
					out.write(newline+'\n')
				else:
					out.write(line)
			elif linenum % 4 == 2:
				length_1 = len(line.strip('\n'))
				out.write(line)
			elif linenum % 4 == 3:
				out.write(line)
			elif linenum % 4 == 0:
				out.write(line)

	else:
		pass
	
	out.close()
	end1.close()
	end2.close()
	singles.close()

if __name__ == "__main__":
		main(sys.argv[1:])