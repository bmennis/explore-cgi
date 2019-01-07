"""Annotations for kaviar with cgi only indels greater than 10bp away from illumina only."""
include: 'const.py'
from collections import defaultdict
import csv

rule ind_split_kaviar_chrom:
    input:  i = DATA + 'interim/kaviar.vcf'
    output: expand(DATA + 'interim/cgi_ind_exp/kaviar/{chrom}.vcf', chrom=list(range(1,23)) + ['X', 'Y', 'M'])
    run:
        outs = {}
        for chrom in list(range(1,23)) + ['X', 'Y', 'M']:
            outs[str(chrom)] = open(DATA + 'interim/cgi_ind_exp/kaviar/' + str(chrom) + '.vcf', 'w')
        with open(input.i) as f:
            for line in f:
                if line[0] == '#':
                    for chrom in outs:
                        print(line.strip(), file=outs[chrom])
                else:
                    chrom = line.split('\t')[0]
                    print(line.strip(), file=outs[chrom])

rule ind_kaviar_vcfanno:
    input:   vcf = DATA + 'interim/cgi_ind_exp/kaviar/{chr}.vcf',
             config = CONFIG + 'kaviar_vcfanno.conf'
    output:  DATA + 'interim/cgi_ind_exp/kaviar_anno/{chr}.vcf'
    threads: 3
    shell:   'vcfanno -p {threads} -base-path {GEMINI_DIR} {input.config} {input.vcf} > {output}'

rule all_ind_kaviar_vcfanno:
    input: expand(DATA + 'interim/cgi_ind_exp/kaviar_anno/{chr}.vcf', chr=list(range(1,23)) + ['X','Y'])

rule ind_parse_kaviar_anno:
    input:  DATA + 'interim/cgi_ind_exp/kaviar_anno/{chr}.vcf'
    output: DATA + 'interim/cgi_ind_exp/kaviar_anno_parsed/{chr}.mat'
    shell:  'python {SCRIPTS}mk_kaviar_matrix.py {input} {output}'

rule ind_mk_kaviar_bed:
    input:  i = DATA + 'interim/cgi_ind_exp/kaviar_anno_parsed/{chr}.mat'
    output:
        cgi = DATA + 'interim/cgi_ind_exp/kaviar_bed/{chr}.cgi.{var_type}.bed',
        both = DATA + 'interim/cgi_ind_exp/kaviar_bed/{chr}.both.{var_type}.bed',
        ill = DATA + 'interim/cgi_ind_exp/kaviar_bed/{chr}.ill.{var_type}.bed',
        cgi_mat = DATA + 'interim/cgi_ind_exp/kaviar_mat/{chr}.cgi.{var_type}.mat',
        both_mat = DATA + 'interim/cgi_ind_exp/kaviar_mat/{chr}.both.{var_type}.mat',
        ill_mat = DATA + 'interim/cgi_ind_exp/kaviar_mat/{chr}.ill.{var_type}.mat'
    run:
        dat = pd.read_csv(input.i, sep='\t')
        dat = dat[dat.var_type==wildcards.var_type]
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

rule all_ind_mk_kaviar_bed:
    input: 
        expand(DATA + 'interim/cgi_ind_exp/kaviar_bed/{chr}.cgi.{var_type}.bed', chr=list(range(1,23)) + ['X','Y'], var_type=('indel')),
        expand(DATA + 'interim/cgi_ind_exp/kaviar_bed/{chr}.both.{var_type}.bed', chr=list(range(1,23)) + ['X','Y'], var_type=('indel')),
        expand(DATA + 'interim/cgi_ind_exp/kaviar_bed/{chr}.ill.{var_type}.bed', chr=list(range(1,23)) + ['X','Y'], var_type=('indel')),
        expand(DATA + 'interim/cgi_ind_exp/kaviar_mat/{chr}.cgi.{var_type}.mat', chr=list(range(1,23)) + ['X','Y'], var_type=('indel')),
        expand(DATA + 'interim/cgi_ind_exp/kaviar_mat/{chr}.both.{var_type}.mat', chr=list(range(1,23)) + ['X','Y'], var_type=('indel')),
        expand(DATA + 'interim/cgi_ind_exp/kaviar_mat/{chr}.ill.{var_type}.mat', chr=list(range(1,23)) + ['X','Y'], var_type=('indel'))

rule ind_mk_kaviar_fasta:
    input:  bed = DATA + 'interim/cgi_ind_exp/kaviar_bed/{chr}.{set}.bed',
            fa = HG19_FA_NOCHR
    output: DATA + 'interim/cgi_ind_exp/kaviar_fa/{chr}.{set}.fa'
    shell:  'bedtools getfasta -fi {input.fa} -bed {input.bed} > {output}'

rule ind_upper_fa:
    input:  i = DATA + 'interim/cgi_ind_exp/kaviar_fa/{aset}.fa'
    output: o = DATA + 'interim/cgi_ind_exp/kaviar_fa_upper/{aset}.fa'
    run:
        with open(input.i) as f, open(output.o, 'w') as fout:
            for line in f:
                if line[0] == '>':
                    print(line.strip(), file=fout)
                else:
                    upper_nucs = line.upper().strip()
                    print(upper_nucs, file=fout)
                    #pur_py = lambda x: 'Y' if x in ('C','T') else 'R'
                    #cg = lambda x: 'S' if x in ('C','G') else 'W'
                    #print(''.join([pur_py(x) for x in upper_nucs]), file=fout)
                    #print(''.join([cg(x) for x in upper_nucs]), file=fout)

rule ind_collapse_kaviar_fasta:
    input:  expand(DATA + 'interim/cgi_ind_exp/kaviar_fa_upper/{chr}.{{aset}}.fa', chr=list(range(1,23)) + ['X', 'Y',])
    output: DATA + 'interim/cgi_ind_exp/kaviar_fa_gz/{aset}.fa.gz'
    shell:  'cat {input} | gzip - > {output}'

rule ind_all_collapse_kav_fa:
    input: expand(DATA + 'interim/cgi_ind_exp/kaviar_fa_gz/{aset}.fa.gz', aset=('cgi.indel','both.indel'))

#rule mk_kaviar_filt_fasta:
#    input: bed = DATA + 'interim/cgi_ind_exp/kaviar_pysster_bed/{chr}.{aset}.bed',
#           fa = HG19_FA_NOCHR
#    output: DATA + 'interim/cgi_ind_exp/kaviar_pysster_fa_ref/{chr}.{aset}.fa'
#    shell: 'bedtools getfasta -fi {input.fa} -bed {input.bed} > {output}'
#
#rule ind_upper_kaviar_filt_fa:
#    input:  i = DATA + 'interim/cgi_ind_exp/kaviar_pysster_fa_ref/{aset}.fa'
#    output: o = DATA + 'interim/cgi_ind_exp/kaviar_pysster_fa_ref_upper/{aset}.fa'
#    run:
#        with open(input.i) as f, open(output.o, 'w') as fout:
#            for line in f:
#                if line[0] == '>':
#                    print(line.strip(), file=fout)
#                else:
#                    upper_nucs = line.upper().strip()
#                    print(upper_nucs, file=fout)
#
#rule ind_collapse_kaviar_filt_fasta:
#    input:  expand(DATA + 'interim/cgi_ind_exp/kaviar_pysster_fa_ref_upper/{chr}.{{aset}}.fa', chr=list(range(1,23)) + ['X', 'Y',])
#    output: DATA + 'interim/cgi_ind_exp/kaviar_pysster_fa_ref_gz/{aset}.fa.gz'
#    shell:  'cat {input} | gzip - > {output}'
#
#rule ind_all_collapse_kav_filt_fa:
#    input: expand(DATA + 'interim/cgi_ind_exp/kaviar_pysster_fa_ref_gz/{aset}.fa.gz', aset=('cgi.indel','both.indel'))

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

rule ind_mk_kaviars:
    input:  #s = expand(DATA + 'interim/cgi_ind_exp/kaviar_anno_parsed/{chr}.full.mat', chr=list(range(1,23)) + ['X', 'Y',]),
            cgi_indel = expand(DATA + 'interim/cgi_ind_exp/kaviar_mat_filt_flanked/{chr}.indel.cgi.filt.flanked.mat', chr=list(range(1,23)) + ['X','Y']),
            both_indel = expand(DATA + 'interim/cgi_ind_exp/kaviar_mat_filt_flanked/{chr}.indel.both.filt.flanked.mat', chr=list(range(1,23)) + ['X','Y']),
            #cgi_snv = expand(DATA + 'interim/cgi_ind_exp/kaviar_mat/{chr}.cgi.snv.mat', chr=list(range(1,23)) + ['X','Y']),
            #both_snv = expand(DATA + 'interim/cgi_ind_exp/kaviar_mat/{chr}.both.snv.mat', chr=list(range(1,23)) + ['X','Y'])
    output: #o = DATA + 'interim/cgi_ind_exp/kaviar.mat',
            cgi_only_indel_out = DATA + 'interim/cgi_ind_exp/cat_mat/cgi_indel.mat',
            both_only_indel_out = DATA + 'interim/cgi_ind_exp/cat_mat/both_indel.mat',
            cgi_indel_out = DATA + 'interim/cgi_ind_exp/kaviar_cgi_indel.mat',
            #both_indel_out = DATA + 'interim/cgi_ind_exp/kaviar_both_indel.mat',
            #cgi_snv_out = DATA + 'interim/cgi_ind_exp/kaviar_cgi_snv.mat',
            #both_snv_out = DATA + 'interim/cgi_ind_exp/kaviar_both_snv.mat'
    run:
        #write_mat(input.s,output.o)
        write_mat(input.cgi_indel, output.cgi_only_indel_out)
        write_mat(input.both_indel, output.both_only_indel_out)
        write_mat(input.cgi_indel + input.both_indel,output.cgi_indel_out)
        #write_mat(input.both_indel,output.both_indel_out)
        #write_mat(input.cgi_snv,output.cgi_snv_out)
        #write_mat(input.both_snv,output.both_snv_out)

rule mk_sample_mat:
    input: DATA + 'interim/cgi_ind_exp/cat_mat/{type}_{var_type}.mat'
    output: 
        samp_mat = DATA + 'interim/cgi_ind_exp/pysster_samp_mat/{type}.{var_type}.sample.mat',
        unsamp_mat = DATA + 'interim/cgi_ind_exp/pysster_unsamp_mat/{type}.{var_type}.unsample.mat'
    shell: 'python {SCRIPTS}sample_mat.py {input} {output.samp_mat} {output.unsamp_mat}'

rule mk_all_sample_mat:
    input: 
        samp_mat = expand(DATA + 'interim/cgi_ind_exp/pysster_samp_mat/{type}.{var_type}.sample.mat', type = ('cgi','both'), var_type = ('indel')),
        unsamp_mat = expand(DATA + 'interim/cgi_ind_exp/pysster_unsamp_mat/{type}.{var_type}.unsample.mat', type = ('cgi','both'), var_type = ('indel'))

rule mk_pos_beds:
    input: DATA + 'interim/cgi_ind_exp/kaviar_mat/{chr}.{s_type}.{var_type}.mat'
    output: DATA + 'interim/cgi_ind_exp/kaviar_bed_pos/{chr}.{s_type}.{var_type}.pos.bed'
    shell: '''awk -F"\t" '{{if (NR!=1) print $4 "\t" $5 "\t" $5+1}}' {input} > {output}'''

rule ind_sort_beds:
    input: DATA + 'interim/cgi_ind_exp/kaviar_bed_pos/{chr}.{s_type}.{var_type}.pos.bed'
    output: DATA + 'interim/cgi_ind_exp/kaviar_bed_sorted/{chr}.{s_type}.{var_type}.sort.bed'
    shell: 'sortBed -i {input} > {output}'

rule ind_mk_closest_beds:
    input: a = DATA + 'interim/cgi_ind_exp/kaviar_bed_sorted/{chr}.cgi.{var_type}.sort.bed',
           b_ill = DATA + 'interim/cgi_ind_exp/kaviar_bed_sorted/{chr}.ill.{var_type}.sort.bed',
           b_both = DATA + 'interim/cgi_ind_exp/kaviar_bed_sorted/{chr}.both.{var_type}.sort.bed'
    output: DATA + 'interim/cgi_ind_exp/kaviar_bed_closest/{chr}.{var_type}.cgi.closest.bed'
    shell: 'bedtools closest -a {input.a} -b {input.b_ill} {input.b_both} -mdb all -d  > {output}'

#rule cat_closest_beds:
#    input: expand(DATA + 'interim/kaviar_bed_closest/{chr}.{{set}}.cgi.closest.bed', chr=list(range(1,23)) + ['X', 'Y',])
#    output: DATA + 'interim/kaviar_bed_closest_cat/{set}.cgi.closest.bed'
#    shell: 'cat {input} > {output}'

#rule all_closest_beds:
#    input: expand(DATA + 'interim/kaviar_bed_closest_cat/{set}.cgi.closest.bed', set=('full.ill.indel', 'full.both.indel'))

rule ind_filter_closest_beds:
    input: DATA + 'interim/cgi_ind_exp/kaviar_bed_closest/{chrom}.{var_type}.cgi.closest.bed'
    output: DATA + 'interim/cgi_ind_exp/kaviar_bed_closest_filtered/{chrom}.{var_type}.cgi.closest.filtered.bed'
    shell: '''awk '{{if (($8 > 0) && ($8 < 10)) print $1 "\t" $2}}' {input} > {output}'''

#rule removed_ind:
#    input: expand(DATA + 'interim/cgi_ind_exp/kaviar_bed_closest/{chrom}.{s_type}.{{var_type}}.cgi.closest.bed', s_type=('ill','both'), chrom=list(range(1,23)) + ['X','Y'])
#    output: DATA + 'interim/cgi_ind_exp/kaviar_bed_closest_removed/{var_type}.cgi.removed.out'
#    shell: '''cat {input} | awk '{{if (($7 > 0) && ($7 < 10)) print $1 "\t" $2}}' | sort -u > {output}'''

#rule all_removed_ind:
#    input: expand(DATA + 'interim/cgi_ind_exp/kaviar_bed_closest_removed/{var_type}.cgi.removed.out', var_type=('indel'))

rule indel_dict:
    input: bed = DATA + 'interim/cgi_ind_exp/kaviar_bed_closest_filtered/{chrom}.{var_type}.cgi.closest.filtered.bed',
           mat = DATA + 'interim/cgi_ind_exp/kaviar_mat/{chrom}.cgi.indel.mat'
    output: DATA + 'interim/cgi_ind_exp/kaviar_mat_filtered/{chrom}.{var_type}.cgi.filtered.mat'
    run:
        cgi_ind = defaultdict(dict)
        with open(input.bed, 'r') as f:
            for line in f:
                pos = line.split('\t')[1]
                pos = pos.strip('\n')
                cgi_ind[wildcards.chrom][pos] = True
        print('100017487' in cgi_ind[wildcards.chrom])
        with open(input.mat, 'r') as f, open(str(output), 'w') as fout:
            reader = csv.DictReader(f, delimiter='\t')
            print('\t'.join(reader.fieldnames), file=fout)
            for row in reader:
                new_pos = row['pos']
                if '100017487' == new_pos:
                  print('found')
                #print(new_pos)
                if new_pos in cgi_ind[wildcards.chrom]:
                    pass
                    #row['closest'] = '1'
                    #print('\t'.join([row[x] for x in reader.fieldnames]), file=fout)
                else:
                    print('\t'.join([row[x] for x in reader.fieldnames]),file=fout)

rule tmp:
  input: DATA + 'interim/cgi_ind_exp/kaviar_mat_filtered/10.indel.cgi.filtered.mat'

rule ind_all_filtered_mat:
    input: expand(DATA + 'interim/cgi_ind_exp/kaviar_mat_filtered/{chrom}.{var_type}.cgi.filtered.mat', chrom=list(range(1,23)) + ['X','Y'], var_type=('indel'))

rule ind_flank_mat:
    input: fa = DATA + 'interim/cgi_ind_exp/kaviar_fa_upper/{chr}.cgi.{var_type}.fa',
           mat = DATA + 'interim/cgi_ind_exp/kaviar_mat_filtered/{chr}.{var_type}.cgi.filtered.mat'
    output: DATA + 'interim/cgi_ind_exp/kaviar_mat_filt_flanked/{chr}.{var_type}.cgi.filt.flanked.mat'
    shell: 'python {SCRIPTS}mk_flanking_seq.py {input.fa} {input.mat} {output}'

rule ind_flank_both:
    input: fa = DATA + 'interim/cgi_ind_exp/kaviar_fa_upper/{chr}.both.{var_type}.fa',
           mat = DATA + 'interim/cgi_ind_exp/kaviar_mat/{chr}.both.{var_type}.mat'
    output: DATA + 'interim/cgi_ind_exp/kaviar_mat_filt_flanked/{chr}.{var_type}.both.filt.flanked.mat'
    shell: 'python {SCRIPTS}mk_flanking_seq.py {input.fa} {input.mat} {output}'

rule ind_mk_samp_add_feat:
    input: expand(DATA + 'interim/cgi_ind_exp/pysster_samp_mat/{type}.{var_type}.sample.mat', type=('cgi','both'), var_type=('indel'))
    shell: 'python {SCRIPTS}mk_pysster_add_feats.py {input}'

rule ind_mk_unsamp_add_feat:
    input: expand(DATA + 'interim/cgi_ind_exp/pysster_unsamp_mat/{type}.{var_type}.unsample.mat', type=('cgi','both'), var_type=('indel'))
    shell: 'python {SCRIPTS}mk_pysster_add_feats.py {input}'

rule ind_mk_filt_beds:
    input: i = DATA + 'interim/cgi_ind_exp/kaviar_mat_filt_flanked/{chr}.{var_type}.{type}.filt.flanked.mat'
    output: o = DATA + 'interim/cgi_ind_exp/kaviar_pysster_bed/{chr}.{type}.{var_type}.bed'
    run:
        dat = pd.read_csv(input.i, sep='\t')
        dat_closest = dat[(dat['closest'] == 0) & (dat['pos_minus'] > 1)]
        dat_closest[['chrom', 'pos_minus', 'pos_plus']].to_csv(output.o, index=False, header=False, sep='\t')

rule ind_mk_all_filt_beds:
    input: expand(DATA + 'interim/cgi_ind_exp/kaviar_pysster_bed/{chr}.{type}.{var_type}.bed', type=('cgi','both'), chr=list(range(1,23)) + ['X','Y'], var_type=('indel'))

rule ind_mk_pysster_sample_fa:
    input: DATA + 'interim/cgi_ind_exp/pysster_samp_mat/{type}.{var_type}.sample.mat'
    output: DATA + 'interim/cgi_ind_exp/pysster_fa/{type}.{var_type}.sample.fa'
    shell: 'python {SCRIPTS}mk_pysster_indel_fa.py {input} {output}'

rule ind_mk_pysster_unsamp_fa:
    input: DATA + 'interim/cgi_ind_exp/pysster_unsamp_mat/{type}.{var_type}.unsample.mat'
    output: DATA + 'interim/cgi_ind_exp/pysster_fa/{type}.{var_type}.unsample.fa'
    shell: 'python {SCRIPTS}mk_pysster_indel_fa.py {input} {output}'

rule ind_all_pysster_fa:
    input: expand(DATA + 'interim/cgi_ind_exp/pysster_fa/{type}.{var_type}.{samp_type}.fa', type=('cgi','both'), var_type=('indel'), samp_type=('sample','unsample'))

rule ind_collapse_pysster_fa:
    input: DATA + 'interim/cgi_ind_exp/pysster_fa/{type}.{var_type}.{samp_type}.fa'
    output: DATA + 'interim/cgi_ind_exp/pysster_fa_gz/{type}.{var_type}.{samp_type}.fa.gz'
    shell: 'gzip {input} > {output}'

rule ind_all_collapse_fa:
    input: 
        fa = expand(DATA + 'interim/cgi_ind_exp/pysster_fa_gz/{type}.{var_type}.{samp_type}.fa.gz', type=('cgi','both'), var_type=('indel'), samp_type=('sample','unsample'))



rule ind_mk_pysster_sample_ref_fa:
    input: DATA + 'interim/cgi_ind_exp/pysster_samp_mat/{type}.{var_type}.sample.mat'
    output: DATA + 'interim/cgi_ind_exp/pysster_ref_fa/{type}.{var_type}.sample.fa'
    shell: 'python {SCRIPTS}mk_pysster_indel_ref_fa.py {input} {output}'

rule ind_mk_pysster_unsamp_ref_fa:
    input: DATA + 'interim/cgi_ind_exp/pysster_unsamp_mat/{type}.{var_type}.unsample.mat'
    output: DATA + 'interim/cgi_ind_exp/pysster_ref_fa/{type}.{var_type}.unsample.fa'
    shell: 'python {SCRIPTS}mk_pysster_indel_ref_fa.py {input} {output}'

rule ind_all_pysster_ref_fa:
    input: expand(DATA + 'interim/cgi_ind_exp/pysster_ref_fa/{type}.{var_type}.{samp_type}.fa', type=('cgi','both'), var_type=('indel'), samp_type=('sample','unsample'))

rule ind_collapse_pysster_ref_fa:
    input: DATA + 'interim/cgi_ind_exp/pysster_ref_fa/{type}.{var_type}.{samp_type}.fa'
    output: DATA + 'interim/cgi_ind_exp/pysster_ref_fa_gz/{type}.{var_type}.{samp_type}.fa.gz'
    shell: 'gzip {input} > {output}'

rule ind_all_collapse_pysster_ref_fa:
    input:
        fa = expand(DATA + 'interim/cgi_ind_exp/pysster_ref_fa_gz/{type}.{var_type}.{samp_type}.fa.gz', type=('cgi','both'), var_type=('indel'), samp_type=('sample','unsample'))



