import argparse, sys
import pandas as pd
sys.path.append('/home/ennisb/me/explore-cgi/src/rules/')
import const

anno_beds = const.ANNO_BEDS + ['ahmad_status', 'indel_length', 'closest']

def main(args):
    files = [args.cgiMatFile, args.bothMatFile]
    for input_file in files:
        i_file = str(input_file)
        with open(input_file, 'r') as inf:
            df = pd.read_csv(inf, sep='\t')
            crit = df.apply(lambda row: row['alt'][0] == row['ref'][0], axis=1)
            df_new = df[crit] # selects True rows
            sp_type = i_file.split('/')[-1].split('.')[0]
            of = i_file.split('/')[0:-2] + ['add_feat'] + [sp_type]
            of = '/'.join(of)
            for flag in anno_beds:
                df1 = df_new[[flag]]
                ofw = of + '__' + flag + '.out'
                with open(ofw, 'w') as outf:
                    df1.to_csv(outf, index=False, header=False)

if __name__ == "__main__":
    desc = 'Mk add feat files for pysster from mat files for individual platforms.'
    parser = argparse.ArgumentParser(description=desc)
    argLs = ('cgiMatFile', 'bothMatFile',)
    for param in argLs:
        parser.add_argument(param)
    args = parser.parse_args()
    main(args) 

