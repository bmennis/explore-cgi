import sys
import csv
import random
#pos 299 in flank seq appears to be indel site

def mk_fa_seq(mat_file, outFile):
    with open(mat_file, 'r') as inf, open(outFile, 'w') as of:
        reader = csv.DictReader(inf, delimiter='\t')
        for row in reader:
            new_fa = row['flanking sequence']
            header = '>' + row['chrom'] + ':' + row['pos_minus'] + '-' +  row['pos_plus'] + '\n'
            print(header + new_fa, file=of)
            ###add back in ind_seq to print statement, and add \n to end of seq

def main(matFile, outFile):
    #reader = mk_df(matFile)
    mk_fa_seq(matFile, outFile)

if __name__ == "__main__":
    matFile = sys.argv[1]
    outFile = sys.argv[2]
    main(matFile, outFile)
