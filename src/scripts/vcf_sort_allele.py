import gzip
import csv
import sys
import os

def checkHeader(outfile):
	"""A function to check if the header is present in outfile and if missing, call the write header function."""

	with open(outfile, 'r') as ofh:
        	header = ofh.readline()
        	if header == '':
                	ofh.close()
                	writeHeader(outfile)

        	else:
                	if header[0] == 'S':
                        	ofh.close()
                	else:
                        	ofh.close()
                        	writeHeader(outfile)

def writeInfo(outfile, source, alleleFrq):
	"""A function to write information to a vcf file."""

	with open(outfile, 'a') as ofh:
		rowWrite = source + '\t' + alleleFrq +  '\n'
		ofh.write(rowWrite)

def writeHeader(outfile):
	"""A function to write the header to the outfile if no header already exists."""
	
	with open(outfile, 'w') as ofh:
		header = "Sources\tAllele Frequency\n"
		ofh.write(header)

def getSource(infile):
	"""A function to extract allele frequency source name from input file name."""
	
	source = os.path.splitext(os.path.basename(infile))[0]
	return source

def parseCsvFile(fileHandle, outfile, source):
	"""A function to parse through the csv file and extract allele frequencies."""
	
	readFile = csv.reader(fileHandle, delimiter='\t')  #read the vcf file using csv module and specifying tab delimitation
	for row in readFile:
		if row[0][0] == '#':    #skip heading lines by checking for #
			next(readFile)
		else:
			info = row[7].split(';')  #split the info section of vcf file by semicolon to read data source
			for x in info:
				if x[0:2] == 'AF':    #check each entry in info section to find allele frequency information
					alleleFrq = x[3:]    #slice the allele frequency cutting out the AF= and grab just the values
					writeInfo(outfile, source, alleleFrq) #write the source and allele frequency to the file entered

				else:
					pass

#Establish files to write variants to after classification/sorting
infile = sys.argv[1]
outfile = sys.argv[2]

source = getSource(infile)
checkHeader(outfile)

if str(sys.argv[1]).endswith(".gz") or str(sys.argv[1]).endswith(".gzip"):
	with gzip.open(infile, 'rt') as gfh:
		parseCsvFile(gfh, outfile, source)

else:
	with open(infile, 'r') as fh:  #open the vcf file in read mode specified
		parseCsvFile(fh, outfile, source)

#		readFile = csv.reader(fh, delimiter='\t')  #read the vcf file using csv module and specifying tab delimitation
	

#	for row in readFile:
		#	if row[0][0] == '#':    #skip heading lines by checking for #
		#		next(readFile)
		#	else:
			#	info = row[7].split(';')  #split the info section of vcf file by semicolon to read data source
			#	for x in info:
				#	if x[0:2] == 'AF':    #check each entry in info section to find allele frequency information
					#	alleleFrq = x[3:]    #slice the allele frequency cutting out the AF= and grab just the values
					#	writeInfo(outfile, source, alleleFrq) #write the source and allele frequency to the file entered

					#else:
						#pass
