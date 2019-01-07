"""Annotations for kaviar with cgi only snvs."""
include: 'const.py'
from collections import defaultdict
import csv

rule snv_mk_kaviar_bed:
    input:  i = DATA + 'interim/cgi_ind_exp/kaviar_anno_parsed/{chr}.mat'
    output:
        cgi = DATA + 'interim/cgi_snv_exp/kaviar_bed/{chr}.cgi.snv.bed',
        both = DATA + 'interim/cgi_snv_exp/kaviar_bed/{chr}.both.snv.bed',
        ill = DATA + 'interim/cgi_snv_exp/kaviar_bed/{chr}.ill.snv.bed',
        cgi_mat = DATA + 'interim/cgi_snv_exp/kaviar_mat/{chr}.cgi.snv.mat',
        both_mat = DATA + 'interim/cgi_snv_exp/kaviar_mat/{chr}.both.snv.mat',
        ill_mat = DATA + 'interim/cgi_snv_exp/kaviar_mat/{chr}.ill.snv.mat'
    run:
        dat = pd.read_csv(input.i, sep='\t')
        dat = dat[dat.var_type=='snv']
        features = [x for x in dat.columns if not x in ('kaviar_status', 'var_type', 'chrom', 'pos', 'ref', 'alt')]
        dat = dat.drop_duplicates(subset=['chrom', 'pos'])
        #dat_ahmad = dat[((dat['ahmad_status'] == 1) | (dat['ahmad_status'] == 0)) ]
        dat_ahmad = dat[dat['ahmad_status'] == 1]
        dat_ahmad.loc[:, 'pos_minus'] = dat_ahmad['pos'] - 300
        dat_ahmad.loc[:, 'pos_plus'] = dat_ahmad['pos'] + 300

        cgi_only_dat = dat_ahmad[(dat_ahmad.kaviar_status=='cgi') & (dat_ahmad.pos_minus>1)]
        both_dat = dat_ahmad[(dat_ahmad.kaviar_status=='both') & (dat_ahmad.pos_minus>1)]
        ill_only_dat = dat_ahmad[(dat_ahmad.kaviar_status=='ill') & (dat_ahmad.pos_minus>1)]
        cgi_only_dat.to_csv(output.cgi_mat, index=False, header=True, sep='\t')
        both_dat.to_csv(output.both_mat, index=False, header=True, sep='\t')
        ill_only_dat.to_csv(output.ill_mat, index=False, header=True, sep='\t')
        cgi_only_dat[['chrom', 'pos_minus', 'pos_plus']].to_csv(output.cgi, index=False, header=False, sep='\t')
        both_dat[['chrom', 'pos_minus', 'pos_plus']].to_csv(output.both, index=False, header=False, sep='\t')
        ill_only_dat[['chrom', 'pos_minus', 'pos_plus']].to_csv(output.ill, index=False, header=False, sep='\t')

rule all_snv_mk_kaviar_bed:
    input:
        expand(DATA + 'interim/cgi_snv_exp/kaviar_bed/{chr}.cgi.snv.bed', chr=list(range(1,23)) + ['X','Y']),
        expand(DATA + 'interim/cgi_snv_exp/kaviar_bed/{chr}.both.snv.bed', chr=list(range(1,23)) + ['X','Y']),
        expand(DATA + 'interim/cgi_snv_exp/kaviar_bed/{chr}.ill.snv.bed', chr=list(range(1,23)) + ['X','Y']),
        expand(DATA + 'interim/cgi_snv_exp/kaviar_mat/{chr}.cgi.snv.mat', chr=list(range(1,23)) + ['X','Y']),
        expand(DATA + 'interim/cgi_snv_exp/kaviar_mat/{chr}.both.snv.mat', chr=list(range(1,23)) + ['X','Y']),
        expand(DATA + 'interim/cgi_snv_exp/kaviar_mat/{chr}.ill.snv.mat', chr=list(range(1,23)) + ['X','Y'])



def write_mat(in_file,out_file):
    with open(list(in_file)[0]) as f:
        header = f.readline().strip()
    with open(out_file, 'w') as fout:
        print(header, file=fout)
        for i in in_file:
            with open(i) as f:
                f.readline()
                for line in f:
                    print(line.strip(), file=fout)
    return out_file

rule snv_mk_kaviars:
    input:  #s = expand(DATA + 'interim/cgi_ind_exp/kaviar_anno_parsed/{chr}.full.mat', chr=list(range(1,23)) + ['X', 'Y',]),
            cgi_snv = expand(DATA + 'interim/cgi_snv_exp/kaviar_mat/{chr}.cgi.snv.mat', chr=list(range(1,23)) + ['X','Y']),
            both_snv = expand(DATA + 'interim/cgi_snv_exp/kaviar_mat/{chr}.both.snv.mat', chr=list(range(1,23)) + ['X','Y']),
            ill_snv = expand(DATA + 'interim/cgi_snv_exp/kaviar_mat/{chr}.ill.snv.mat', chr=list(range(1,23)) + ['X','Y'])
    output: #o = DATA + 'interim/cgi_ind_exp/kaviar.mat',
            cgi_snv_out = DATA + 'interim/cgi_snv_exp/cat_mat/kaviar_cgi_snv.mat',
            both_snv_out = DATA + 'interim/cgi_snv_exp/cat_mat/kaviar_both_snv.mat',
            ill_snv_out = DATA + 'interim/cgi_snv_exp/cat_mat/kaviar_ill_snv.mat'
    run:
        #write_mat(input.s,output.o)
        write_mat(input.cgi_snv,output.cgi_snv_out)
        write_mat(input.both_snv,output.both_snv_out)
        write_mat(input.ill_snv,output.ill_snv_out)
