"""Annotations for kaviar with cgi only indels greater than 10bp away from illumina only."""

from collections import defaultdict

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
            cgi_indel = expand(DATA + 'interim/cgi_ind_exp/kaviar_mat_filtered/{chr}.indel.cgi.filtered.mat', chr=list(range(1,23)) + ['X','Y']),
            both_indel = expand(DATA + 'interim/cgi_ind_exp/kaviar_mat/{chr}.both.indel.mat', chr=list(range(1,23)) + ['X','Y']),
            #cgi_snv = expand(DATA + 'interim/cgi_ind_exp/kaviar_mat/{chr}.cgi.snv.mat', chr=list(range(1,23)) + ['X','Y']),
            #both_snv = expand(DATA + 'interim/cgi_ind_exp/kaviar_mat/{chr}.both.snv.mat', chr=list(range(1,23)) + ['X','Y'])
    output: #o = DATA + 'interim/cgi_ind_exp/kaviar.mat',
            cgi_indel_out = DATA + 'interim/cgi_ind_exp/kaviar_cgi_indel.mat',
            #both_indel_out = DATA + 'interim/cgi_ind_exp/kaviar_both_indel.mat',
            #cgi_snv_out = DATA + 'interim/cgi_ind_exp/kaviar_cgi_snv.mat',
            #both_snv_out = DATA + 'interim/cgi_ind_exp/kaviar_both_snv.mat'
    run:
        #write_mat(input.s,output.o)
        write_mat(input.cgi_indel + input.both_indel,output.cgi_indel_out)
        #write_mat(input.both_indel,output.both_indel_out)
        #write_mat(input.cgi_snv,output.cgi_snv_out)
        #write_mat(input.both_snv,output.both_snv_out)

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
           b = DATA + 'interim/cgi_ind_exp/kaviar_bed_sorted/{chr}.{s_type}.{var_type}.sort.bed'
    output: DATA + 'interim/cgi_ind_exp/kaviar_bed_closest/{chr}.{s_type}.{var_type}.cgi.closest.bed'
    shell: 'bedtools closest -a {input.a} -b {input.b} -d  > {output}'

#rule cat_closest_beds:
#    input: expand(DATA + 'interim/kaviar_bed_closest/{chr}.{{set}}.cgi.closest.bed', chr=list(range(1,23)) + ['X', 'Y',])
#    output: DATA + 'interim/kaviar_bed_closest_cat/{set}.cgi.closest.bed'
#    shell: 'cat {input} > {output}'

#rule all_closest_beds:
#    input: expand(DATA + 'interim/kaviar_bed_closest_cat/{set}.cgi.closest.bed', set=('full.ill.indel', 'full.both.indel'))

rule ind_filter_closest_beds:
    input: expand(DATA + 'interim/cgi_ind_exp/kaviar_bed_closest/{{chrom}}.{s_type}.{{var_type}}.cgi.closest.bed', s_type=('ill','both'))
    output: DATA + 'interim/cgi_ind_exp/kaviar_bed_closest_filtered/{chrom}.{var_type}.cgi.closest.filtered.out'
    shell: '''cat {input} | awk '{{if (($7 == 0) || ($7 >= 10)) print $1 "\t" $2}}' | sort -u > {output}'''

rule indel_dict:
    input: bed = DATA + 'interim/cgi_ind_exp/kaviar_bed_closest_filtered/{chrom}.{var_type}.cgi.closest.filtered.out',
           mat = DATA + 'interim/cgi_ind_exp/kaviar_mat/{chrom}.cgi.indel.mat'
    output: DATA + 'interim/cgi_ind_exp/kaviar_mat_filtered/{chrom}.{var_type}.cgi.filtered.mat'
    run:
        cgi_ind = defaultdict(dict)
        for chrom in list(range(1,23)) + ['X','Y']:
            with open(DATA + 'interim/cgi_ind_exp/kaviar_bed_closest_filtered/' + str(wildcards.chrom) + '.' + str(wildcards.var_type) + '.cgi.closest.filtered.out', 'r') as f:
                for line in f:
                    pos = line.split('\t')[1]
                    pos = pos.strip('\n')
                    cgi_ind[chrom][pos] = True

        for chrom in list(range(1,23)) + ['X','Y']:
            with open(DATA + 'interim/cgi_ind_exp/kaviar_mat/' + str(wildcards.chrom) + '.cgi.indel.mat', 'r') as f, open(str(output), 'w') as o:
                head = f.readline()
                o.write(head)
                for line in f:
                    new_pos = line.split('\t')[4]
                    if new_pos in cgi_ind[str(chrom)]:
                        o.write(line)
                    else:
                        pass

rule ind_all_filtered_mat:
    input: expand(DATA + 'interim/cgi_ind_exp/kaviar_mat_filtered/{chrom}.{var_type}.cgi.filtered.mat', chrom=list(range(1,23)) + ['X','Y'], var_type=('indel'))
