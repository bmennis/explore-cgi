import argparse, sys
import pandas as pd
sys.path.append('/home/ennisb/me/explore-cgi/src/rules/')
import const
import re

anno_beds = const.ANNO_BEDS + ['ahmad_status', 'indel_length', 'closest']

def get_platform(info):
    ###Function to get platform type from name of input file###
    if 'cgi' in info:
        platform = 'cgi'
    elif 'both' or ('cgi' and 'ill') in info:
        platform = 'both'
    elif 'ill' in info:
        platform = 'ill'
    return platform

def get_var_type(info):
    ###Function to get variant type of matrix from name of input file###
    if 'indel' in info:
        var_type = 'indel'
    elif 'snv' in info:
        var_type = 'snv'
    return var_type

def get_crit(df):
    ###Function to get critical variants from indels filtering out those which the first base do not match###
    crit = df.apply(lambda row: row['alt'][0] == row['ref'][0], axis=1)
    df_new = df[crit] # selects True rows

    return df_new

def get_samp_type(info):
    ###Function to get sample or unsample type from name of input file###
    if 'unsample' in info:
        samp_type = 'unsamp'
    elif 'sample' in info:
        samp_type = 'samp'

    return samp_type

def main(args):
    files = list(args.input)
    for input_file in files:
        with open(input_file, 'r') as inf:
            df = pd.read_csv(inf, sep='\t')
            info = str(input_file).split('/')[-1]
            info = re.split('[_.]', info)
            print(info)
            platform = get_platform(info)
            var_type = get_var_type(info)
            samp_type = get_samp_type(info)
            if var_type == 'indel':
                df_new = get_crit(df)
            else:
                df_new = df
            of = str(input_file).split('/')[0:-2] + ['add_feats']
            of = '/'.join(of)
            of = of + '/' + platform + '.' + var_type + '.' + samp_type + '.'
            for flag in anno_beds:
                df1 = df_new[[flag]]
                ofw = of + flag + '.out'
                with open(ofw, 'w') as outf:
                    df1.to_csv(outf, index=False, header=False)

if __name__ == "__main__":
    desc = 'Mk add feat files for pysster from mat files for individual platforms.'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('input', nargs = '+')
    args = parser.parse_args()
    main(args) 

