#Sort the line counts of files generated from pipeline and write to csv file.
import os
import sys
import csv

def checkOutfile(outfile):
	"""A function to check if the header is present in outfile and if missing, call the write header function."""

	if os.path.exists(outfile):
		pass

	else:
		writeHeader(outfile)

def writeHeader(outfile):
	"""A function to write the header to the outfile if no header already exists."""

	with open(outfile, 'w') as ofh:
		header = "Subset\tBed Region Name\tCount\n"
		ofh.write(header)

def getSubset(row):
	"""A function to extract source name from input file name."""
	
	s = ' ' 
	subset_group = os.path.basename(row[1]).split('_')
	if len(subset_group) < 2:
		subset = subset_group[0]
	else:
		if subset_group[1] == 'illumina' or subset_group[1] == 'not':
			subset = s.join(subset_group[0:3])
		else:
			subset = subset_group[0]
	return subset

def getBedRegionName(row):
	"""A function to extract bed region name from line information."""
	
	region = os.path.basename(row[1]).split('_')
	if region[1] == 'illumina' or region[1] == 'not':
		bed_region_name = region[3]
	else:
		bed_region_name = region[1]

	return bed_region_name

def writeInfo(outfile, subset, bed_region_name, line_count):
	"""A function to write information to a csv file."""

	with open(outfile, 'a') as ofh:
		rowWrite = subset + '\t' + bed_region_name + '\t' + line_count +  '\n'
		ofh.write(rowWrite)

def parseFile(fileHandle, outfile):
	"""A function to parse through the file and sort information."""
	readFile = csv.reader(fileHandle, delimiter=' ')
	for row in readFile:
		new_row = list(filter(None, row))
		if new_row[1] == 'total':    #skip lines with total line counts 
			pass
		else:
			line_count = new_row[0]  #establish line count for file
			subset = getSubset(new_row)
			if  new_row[1].endswith('.vcf'):
				bed_region_name = 'NA'
				writeInfo(outfile, subset, bed_region_name, line_count)

			elif new_row[1].endswith('difficult.intr'):
				bed_region_name = getBedRegionName(new_row) + ' Difficult'
				writeInfo(outfile, subset, bed_region_name, line_count)

			else:
				bed_region_name = getBedRegionName(new_row)
				writeInfo(outfile, subset, bed_region_name, line_count)

infile = sys.argv[1]
outfile = sys.argv[2]

checkOutfile(outfile)

with open(infile, 'r') as fh:
	parseFile(fh, outfile)
