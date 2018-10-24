import sys
import os
from collections import defaultdict

cgi_pos_ill_clst = defaultdict(dict)
cgi_pos_both_clst = defaultdict(dict)

def mk_dict(infile, dict_input):
    """Function to write position and distance to a dictionary"""
    
    dict_input = defaultdict(dict)
    with open(infile, 'r') as f:
        for line in f:
            pos = line.split('\t')[1]
            dist = int(line.split('\t')[6])
            dict_input[pos]=dist
    return dict_input
    
with open(input.mat, 'r') as f, open(str(output), 'w') as fout:
        reader = csv.DictReader(f, delimiter='\t')
        print('\t'.join(reader.fieldnames), file=fout)
        for row in reader:
            new_pos = row['pos']
            if '100017487' == new_pos:
              print('found')
            #print(new_pos)
            if new_pos in cgi_ind[wildcards.chrom]:
                print('\t'.join([row[x] for x in reader.fieldnames]), file=fout)
            else:
                row['closest'] = '1'
                print('\t'.join([row[x] for x in reader.fieldnames]),file=fout)

def compare_dict(ill_dict, both_dict, outfile)
    """Compare the dictionaries created and pass information to write to output file"""
    
    int_list= []
    with open(outfile, 'w') as of:
        for pos in ill_dict.keys():
            if pos in both_dict.keys():
                int_list.append(pos)        
            else: 
                if pos[dist] > 0 and pos[dist] < 10:
                    write
                else:
                    pass
        for pos in int_list:
            
