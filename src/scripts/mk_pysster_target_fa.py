import sys
import csv
import random
#pos 299 in flank seq appears to be indel site

#def mk_df(mat_file):
#    with open(mat_file, 'r') as inf:
#        reader = csv.DictReader(inf, delimiter='\t')
#        print(reader)
#    return reader

###TA + 'interim/cgi_ind_exp/pysster_fa/{type}.{var_type}.{samp_type}.af.{af}.fa'

def mk_fa_seq(mat_file, outFile):
    with open(mat_file, 'r') as inf, open(outFile, 'w') as of:
        reader = csv.DictReader(inf, delimiter='\t')
        for row in reader:
            if row['alt'][0] == row['ref'][0]:
                if int(row['indel_length']) > 0:
                    new_fa = row['flanking sequence'][270:300] + row['alt'][1:] + row['flanking sequence'][300:330]
                    ind_seq = ('X' * 30) + ('I' * int(row['indel_length'])) + ('X' * 30)
            
                    if int(row['indel_length']) % 2 == 0:
                        cut = int(int(row['indel_length']) / 2)
                        row['seq_fa'] = new_fa[cut:-cut]
                        row['ind_seq'] = ind_seq[cut:-cut]

                    else:
                        cut_l = int(int(row['indel_length']) / 2)
                        cut_r = int((int(row['indel_length']) / 2) + 1)
                        row['seq_fa'] = new_fa[cut_l:-cut_r]
                        row['ind_seq'] = ind_seq[cut_l:-cut_r]

                else:
                    row['seq_fa'] = row['flanking sequence'][270:330]
                    row['ind_seq'] = ('X' * 30) + ('D' * abs(int(row['indel_length']))) + ('X' * (60-30-abs(int(row['indel_length']))))
                    if len(row['ind_seq']) > 60:
                        row['ind_seq'] = row['ind_seq'][0:60]
        
                header = '>' + row['chrom'] + ':' + row['pos_minus'] + '-' +  row['pos_plus'] + '\n'
                seq = row['seq_fa'] + '\n'
                ind_seq = row['ind_seq']
                print(header + seq + ind_seq, file=of)
                ###add back in ind_seq to print statement, and add \n to end of seq

            else:
                pass

def main(matFile, outFile):
    #reader = mk_df(matFile)
    mk_fa_seq(matFile, outFile)

if __name__ == "__main__":
    matFile = sys.argv[1]
    outFile = sys.argv[2]
    main(matFile, outFile)
