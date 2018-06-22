import csv
import sys

def write_output(new_line, output_file):
    """A function to write information to output vcf file."""
	
    with open(output_file, 'a') as ofh:
        ofh.write(new_line)

hdr = """
##fileformat=VCFv4.1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n"""

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file, 'r') as ifh:
    readFile = csv.reader(ifh, delimiter='\t')
    write_output(hdr, output_file)
    for row in readFile:            
        new_line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (row[0], row[1], '.', row[2], row[3], '.', '.', '.')
        write_output(new_line, output_file)
