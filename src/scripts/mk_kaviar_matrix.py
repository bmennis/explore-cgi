"""Make matrix from vcf file"""
import argparse, sys
sys.path.append('/home/ennisb/me/explore-cgi/src/rules/')
import const

ms_head = ['ms1', 'ms2', 'ms3', 'ms4', 'ms5', 'ms6', 'ms7', 'ms8', 'ms9', 'ms10', 'ms11', 'ms12']

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

ANNO_BEDS = const.ANNO_BEDS #['lowmappabilityall', 'notinlowmappabilityall', 'siren_similarRegions_dist1']
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

def mk_ahmad_status(line):
    row = line.strip().split('\t')
    if 'ahmad' in row[7]:
        ahmad_status = list(set(row[7].split('ahmad_region=')[1].split(';')[0].split(',')))
        if len(ahmad_status) > 1:
            return ['-1']

        ahmad_status = ahmad_status[0]
        if ahmad_status == 'good':
            return ['1']
        elif ahmad_status == 'poor':
            return ['0']
        else:
            print(ahmad_status)
            i = 1/0
    return ['-1']

def mk_af(line):
    row = line.strip().split('\t')
    af = row[7].split('AF=')[1].split(';')[0]
    return af

def mk_indel_length(line):
    ref, alt = line.split('\t')[3:5]
    ind_len = str(len(alt) - len(ref))
    return [ind_len]

def mk_microsat_info(line):
    row = line.strip().split('\t')
    if 'microsat_info' in row[7]:
        microsat = row[7].split('microsat_info=')[1].split(';')[0]
        if ',' in microsat:
            ms_info = microsat.split(',')[0].split('_') 
        else:
            ms_info = microsat.strip().split('_')
    else:
        ms_info = ['MISS'] * 12
    return ms_info

def main(args):
    with open(args.vcfFile) as f, open(args.outFile, 'w') as fout:
        header = ['af', 'kaviar_status', 'var_type', 'chrom', 'pos', 'ref', 'alt'] + ANNO_BEDS + ['ahmad_status', 'indel_length', 'closest'] + ms_head
        print('\t'.join(header), file=fout)
        for line in f:
            if line[0] != '#':
                kaviar_status = get_kaviar_status(line)
                ahmad_status = mk_ahmad_status(line)
                var_type = mk_var_type(line)
                af = mk_af(line)
                ind_len = mk_indel_length(line)
                ms_info = mk_microsat_info(line)
                ls = [af, kaviar_status, var_type] + mk_line_info(line) + ahmad_status + ind_len + ['0'] + ms_info
                print('\t'.join(ls), file=fout)

if __name__ == "__main__":
    desc = 'Mk matrix from vcf.'
    parser = argparse.ArgumentParser(description=desc)
    argLs = ('vcfFile', 'outFile',)
    for param in argLs:
        parser.add_argument(param)
    args = parser.parse_args()
    main(args)
