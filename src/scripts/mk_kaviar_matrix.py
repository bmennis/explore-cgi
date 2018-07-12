"""Make matrix from vcf file"""
import argparse

def get_kaviar_status(line):
    row = line.strip().split('\t')
    if 'DS=' in row[7]:
        info = row[7].split('DS=')[1].split(';')[0]  #split the info section of vcf file by semicolon to read data source
        sources = info.split('|')    #split the data source info by pipe and sort based on following conditions
        cgi_ls = [x for x in sources if x[0:2] in ('NA','GS') or 'Wellderly' in x or x in ('ISB_founders-Nge3','Inova_CGI_founders-Nge3')]
        ill_ls = [x for x in sources if not x[0:2] in ('NA','GS') and not 'Wellderly' in x and not x in ('ISB_founders-Nge3','Inova_CGI_founders-Nge3')]
        if cgi_ls and ill_ls:
            return 'both'
        elif cgi_ls and not ill_ls:
            return 'cgi'
        elif ill_ls and not cgi_ls:
            return 'ill'
    return 'none'

ANNO_BEDS = ['lowmappabilityall', 'notinlowmappabilityall', 'siren_similarRegions_dist1']
def mk_flags(info_ls):
    flags = []
    for flag in ANNO_BEDS:
        if flag in info_ls:
            flags.append('1')
        else:
            flags.append('0')
    return flags

def mk_line_info(line):
    sp = line.strip().split('\t')
    return sp[0:2] + sp[3:5] + mk_flags(sp[7].split(';'))

def mk_var_type(line):
    ref, alt = line.split('\t')[3:5]
    if len(ref) == 1 and len(alt) == 1:
        return 'snv'
    if len(ref) == len(alt) and ref != '.' and alt != '.':
        return 'subs'
    return 'indel'

def main(args):
    with open(args.vcfFile) as f, open(args.outFile, 'w') as fout:
        header = ['kaviar_status', 'var_type', 'chrom', 'pos', 'ref', 'alt'] + ANNO_BEDS
        print('\t'.join(header), file=fout)
        for line in f:
            if line[0] != '#':
                kaviar_status = get_kaviar_status(line)
                var_type = mk_var_type(line)
                ls = [kaviar_status, var_type] + mk_line_info(line)
                print('\t'.join(ls), file=fout)

if __name__ == "__main__":
    desc = 'Mk matrix from vcf.'
    parser = argparse.ArgumentParser(description=desc)
    argLs = ('vcfFile', 'outFile',)
    for param in argLs:
        parser.add_argument(param)
    args = parser.parse_args()
    main(args)
