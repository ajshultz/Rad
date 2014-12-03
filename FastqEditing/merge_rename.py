#!/usr/bin/env python

import sys, getopt

#Note that this script is for taking the output of process_radtags from Stacks, and formatting paired-end reads from double-digest RADseq to be used for de novo assembly.  Input is read 1 and read 2 for an individual, and possibly one or two extra files with unpaired read 1 or read 2 data. This script will concatenate these files, including an option to re-add the read information to the illumina position info (_1 or _2).  


def main(argv):
	try:
		opts,args = getopt.getopt(argv,'h1:2:s:o:n:u:',)
	except getopt.GetoptError:
		print "merge_rename.py -1 <End 1 Fastq Input File> -2 <End 2 Fastq Input File> -s <Unpaired Fastq Input File end 1> -u <Unpaired Fastq Input File end 2> -o <Output File> -n <T or F, Add _1 or _2 back to read names>"
		sys.exit(2)
			
	inputfile1 = ''
	inputfile2 = ''
	inputfilesingles1 = ''
	inputfilesingles2 = ''
	outputfile = ''
	rename = 'F'

	for opt, arg in opts:
		if opt == "-h":
			print "merge_rename.py -1 <End 1 Fastq Input File> -2 <End 2 Fastq Input File> -s <Unpaired Fastq Input File end 1> -u <Unpaired Fastq Input File end 2> -o <Output File> -n <T or F, Add _1 or _2 back to read names>"
			sys.exit()
		elif opt == "-1":
			inputfile1 = arg
		elif opt == "-2":
			inputfile2 = arg
		elif opt == "-s":
			inputfilesingles1 = arg
		elif opt == "-u":
			inputfilesingles2 = arg
		elif opt == "-o":
			outputfile = arg
		elif opt == "-n":
			rename = arg

		

	end1 = open(inputfile1,"r")
	end2 = open(inputfile2,"r")
	if inputfilesingles1 != '':
		singles1 = open(inputfilesingles1,"r")
	else:
		pass
	if inputfilesingles2 != '':
		singles2 = open(inputfilesingles2,"r")
	else:
		pass
	out = open(outputfile,"w")

	linenum = 0
	
	#Iterate through end 1 file, add all lines to the outfile, and add the read information to the Illumina position info if that option is selected.
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

	#Iterate through end 2 file, add all lines to the outfile, and add the read information to the Illumina position info if that option is selected.
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
			out.write(line)
		if linenum % 4 == 3:
			out.write(line)
		if linenum % 4 == 0:
			out.write(line)
			
			
	#Iterate through end 1 singles file if present, add all lines to the outfile, and add the read information to the Illumina position info if that option is selected.
	
	if inputfilesingles1 != '':
	
		linenum = 0
	
		for line in singles1:
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
		
		singles1.close()

	else:
		pass


	#Iterate through end 2 singles file if present, add all lines to the outfile, and add the read information to the Illumina position info if that option is selected.
	
	if inputfilesingles2 != '':
		
		linenum = 0
	
		for line in singles2:
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
			elif linenum % 4 == 2:
				length_1 = len(line.strip('\n'))
				out.write(line)
			elif linenum % 4 == 3:
				out.write(line)
			elif linenum % 4 == 0:
				out.write(line)
		
		singles2.close()	
				
	else:
		pass
	
	out.close()
	end1.close()
	end2.close()
	
if __name__ == "__main__":
		main(sys.argv[1:])