"""Annotations for TARGET."""
include: 'const.py'
FILE_LIST = []

rule filterTargetForIll:
    input: DATA + 'raw/TARGET.mat'
    output: DATA + 'interim/target_exp/target_wxs_filter/target_mat_sorted_wxsCut_{wxsCut}'
    shell: '''cat {input} | awk -F '\\t' 'FNR > 1 {{if($27 > '{wildcards.wxsCut}') {{print $0}}}}' | cut -f 3,4,5,6,7,8,9,10,24,25,26 | sort -u > {output}'''

rule filterTargetNotIll:
    input: DATA + 'raw/TARGET.mat'
    output: DATA + 'interim/target_exp/target_wxs_filter/target_mat_sorted_wxsCut_none'
    shell: '''cat {input} | awk -F '\\t' '($7=="0") && ($8=="0") && ($9=="0") && ($10=="0") && ($24=="0") && ($25=="0") && ($26=="0") {{print $0}}' | cut -f 3,4,5,6,7,8,9,10,24,25,26 | sort -u > {output}'''

rule filterTargetAll:
    input: 
        cgi = DATA + 'interim/target_exp/target_wxs_filter/target_mat_sorted_wxsCut_none',
        ill = expand(DATA + 'interim/target_exp/target_wxs_filter/target_mat_sorted_wxsCut_{wxsCut}', wxsCut = ('10','15','20'))

rule filter_target_ill_dec_tree:
    input: DATA + 'raw/TARGET.mat'
    output: DATA + 'interim/target_exp/target_dec_tree_filter/target_dec_tree_ill_wxsCut_{wxsCut}.mat'
    shell: '''cat {input} | awk -F '\\t' '{{if($27 > '{wildcards.wxsCut}') {{print $0}}}}' | sort -u > {output}'''

rule filter_target_cgi_dec_tree:
    input: DATA + 'raw/TARGET.mat'
    output: DATA + 'interim/target_exp/target_dec_tree_filter/target_dec_tree_cgi_wxsCut_none.mat'
    shell: '''cat {input} | awk -F '\\t' '{{print $0}} ($7=="0") && ($8=="0") && ($9=="0") && ($10=="0") && ($24=="0") && ($25=="0") && ($26=="0") {{print $0}}' | sort -u > {output}'''

rule filter_target_dec_tree_all:
    input:
        cgi = DATA + 'interim/target_exp/target_dec_tree_filter/target_dec_tree_cgi_wxsCut_none.mat',
        ill = expand(DATA + 'interim/target_exp/target_dec_tree_filter/target_dec_tree_ill_wxsCut_{wxsCut}.mat', wxsCut = ('10', '15', '20'))

rule prepTargetMatchMismatchFiles:
    input: DATA + 'interim/target_exp/target_wxs_filter/{subset}'
    output: DATA + 'interim/target_exp/vcf_prep/{subset}_vcf_prep'
    shell:  'cut -f 1,2,3,4 {input} | sort -k 1,1 -k 2,2n > {output}'

rule prepAllTargetMatchMismatch:
    input: expand(DATA + 'interim/target_exp/vcf_prep/{subset}_vcf_prep', subset=('target_mat_sorted_wxsCut_20','target_mat_sorted_wxsCut_none'))

rule createTargetMatchMismatchVcfFile:
    input: DATA + 'interim/target_exp/vcf_prep/target_mat_sorted_{subset}_vcf_prep'
    output: DATA + 'interim/target_exp/vcf/{subset}.vcf'
    shell: 'python {SCRIPTS}vcf_create.py {input} {output}'

rule createAllMatchMismatchVcfFiles:
    input: expand(DATA + 'interim/target_exp/vcf/{subset}.vcf', subset=('wxsCut_20','wxsCut_none')) 

rule target_vcfanno:
    input:   vcf = DATA + 'interim/target_exp/vcf/{subset}.vcf',
             config = CONFIG + 'kaviar_vcfanno.conf'
    output:  DATA + 'interim/target_exp/vcf_anno/{subset}.vcf'
    threads: 3
    shell:   'vcfanno -p {threads} -base-path {GEMINI_DIR} {input.config} {input.vcf} > {output}'

rule all_target_vcfanno:
    input: expand(DATA + 'interim/target_exp/vcf_anno/{subset}.vcf', subset=('wxsCut_20','wxsCut_none'))

rule target_parse_anno:
    input:  DATA + 'interim/target_exp/vcf_anno/{subset}.vcf'
    output: DATA + 'interim/target_exp/vcf_anno_parsed/{subset}.mat'
    shell:  'python {SCRIPTS}mk_target_matrix.py {input} {output}'

rule target_switch_cgi:
    input: DATA + 'interim/target_exp/vcf_anno_parsed/wxsCut_none.mat'
    output: DATA + 'interim/target_exp/vcf_anno_parsed/cgi.mat'
    shell: 'mv {input} {output}'

rule target_switch_both:
    input: DATA + 'interim/target_exp/vcf_anno_parsed/wxsCut_20.mat'
    output: DATA + 'interim/target_exp/vcf_anno_parsed/both.mat'
    shell: 'mv {input} {output}'

rule all_target_switch:
    input: expand(DATA + 'interim/target_exp/vcf_anno_parsed/{platform}.mat', platform = ('cgi','both'))

rule mk_target_bed:
    input:  i = DATA + 'interim/target_exp/vcf_anno_parsed/{platform}.mat'
    output:
        bed = DATA + 'interim/target_exp/target_bed/{platform}.{var_type}.bed',
        mat = DATA + 'interim/target_exp/target_mat/{platform}.{var_type}.mat',
    run:
        dat = pd.read_csv(input.i, sep='\t')
        dat = dat[dat.var_type==wildcards.var_type]
        features = [x for x in dat.columns if not x in ('platform_status', 'var_type', 'chrom', 'pos', 'ref', 'alt')]
        dat = dat.drop_duplicates(subset=['chrom', 'pos'])
        #dat_ahmad = dat[dat['ahmad_status'] == 1]
        dat.loc[:, 'pos_minus'] = dat['pos'] - 300
        dat.loc[:, 'pos_plus'] = dat['pos'] + 300

        platform_only_dat = dat[(dat.platform_status==wildcards.platform) & (dat.pos_minus>1)]
        platform_only_dat.to_csv(output.mat, index=False, header=True, sep='\t')
        platform_only_dat[['chrom', 'pos_minus', 'pos_plus']].to_csv(output.bed, index=False, header=False, sep='\t')

rule all_mk_target_bed:
    input:
        bed = expand(DATA + 'interim/target_exp/target_bed/{platform}.{var_type}.bed', var_type=('indel', 'snv'), platform=('cgi','both')),
        mat = expand(DATA + 'interim/target_exp/target_mat/{platform}.{var_type}.mat', var_type=('indel', 'snv'), platform=('cgi','both'))

#rule mk_target_cgi_add_feat:
#    input: expand(DATA + 'interim/target_exp/target_mat/cgi.{var_type}.mat', var_type=('indel','snv'))
#    shell: 'python {SCRIPTS}mk_add_target_feats.py {input}'

#rule mk_target_both_add_feat:
#    input: expand(DATA + 'interim/target_exp/target_mat/both.{var_type}.mat', var_type=('indel','snv'))
#    shell: 'python {SCRIPTS}mk_add_target_feats.py {input}'

rule mk_target_fasta:
    input:  bed = DATA + 'interim/target_exp/target_bed/{platform}.{var_type}.bed',
            fa = HG19_FA_NOCHR
    output: DATA + 'interim/target_exp/target_fa/{platform}.{var_type}.fa'
    shell:  'bedtools getfasta -fi {input.fa} -bed {input.bed} > {output}'

rule target_upper_fa:
    input:  i = DATA + 'interim/target_exp/target_fa/{aset}.fa'
    output: o = DATA + 'interim/target_exp/target_fa_upper/{aset}.fa'
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

rule collapse_target_fasta:
    input:  DATA + 'interim/target_exp/target_fa_upper/{aset}.fa'
    output: DATA + 'interim/target_exp/target_fa_gz/{aset}.fa.gz'
    shell:  'gzip {input} - > {output}'

rule all_collapse_target_fa:
    input: expand(DATA + 'interim/target_exp/target_fa_gz/{aset}.fa.gz', aset=('cgi.indel','both.indel','cgi.snv','both.snv'))

rule target_flank_mat:
    input: fa = DATA + 'interim/target_exp/target_fa_upper/{platform}.{var_type}.fa',
           mat = DATA + 'interim/target_exp/target_mat/{platform}.{var_type}.mat'
    output: DATA + 'interim/target_exp/target_flanked_mat/{platform}.{var_type}.flanked.mat'
    shell: 'python {SCRIPTS}mk_flanking_seq.py {input.fa} {input.mat} {output}'

rule all_target_flank_mat:
    input: expand(DATA + 'interim/target_exp/target_flanked_mat/{platform}.{var_type}.flanked.mat', platform=('cgi','both'), var_type=('indel','snv'))

rule sample_target_flank_mat:
    input: DATA + 'interim/target_exp/target_flanked_mat/{platform}.{var_type}.flanked.mat'
    output:
        samp_mat = DATA + 'interim/target_exp/target_flanked_mat_sample/{platform}.{var_type}.flanked.sample.mat',
        unsamp_mat = DATA + 'interim/target_exp/target_flanked_mat_unsample/{platform}.{var_type}.flanked.unsample.mat'
    shell: 'python {SCRIPTS}sample_mat.py {input} {output.samp_mat} {output.unsamp_mat}'

rule mk_target_samp_add_feat:
    input: expand(DATA + 'interim/target_exp/target_flanked_mat_sample/{platform}.{var_type}.flanked.sample.mat', platform=('cgi','both'), var_type=('indel'))
    shell: 'python {SCRIPTS}mk_add_target_feats.py {input}'

rule mk_target_unsamp_add_feat:
    input: expand(DATA + 'interim/target_exp/target_flanked_mat_unsample/{platform}.{var_type}.flanked.unsample.mat', platform=('cgi','both'), var_type=('indel'))
    shell: 'python {SCRIPTS}mk_add_target_feats.py {input}'

rule mk_target_pysster_fa:
    input: DATA + 'interim/target_exp/target_flanked_mat_{samp_type}/{type}.{var_type}.flanked.{samp_type}.mat'
    output: DATA + 'interim/target_exp/target_pysster_fa/{type}.{var_type}.{samp_type}.fa'
    shell: 'python {SCRIPTS}mk_pysster_target_fa.py {input} {output}'

rule all_target_pysster_fa:
    input: expand(DATA + 'interim/target_exp/target_pysster_fa/{type}.{var_type}.{samp_type}.fa', type=('cgi','both'), var_type=('indel'), samp_type = ('sample','unsample'))

rule collapse_target_pysster_fa:
    input: DATA + 'interim/target_exp/target_pysster_fa/{type}.{var_type}.{samp_type}.fa'
    output: DATA + 'interim/target_exp/target_pysster_fa_gz/{type}.{var_type}.{samp_type}.fa.gz'
    shell: 'gzip {input} > {output}'

rule all_collapse_target_pysster_fa:
    input:
        fa = expand(DATA + 'interim/target_exp/target_pysster_fa_gz/{type}.{var_type}.{samp_type}.fa.gz', type=('cgi','both'), var_type=('indel'), samp_type = ('sample','unsample'))

rule target_filter_kav_beds:
    input: DATA + 'interim/cgi_ind_exp/kaviar_bed_closest/{chrom}.{var_type}.cgi.closest.bed'
    output: DATA + 'interim/target_exp/kaviar_bed_filtered/{chrom}.{var_type}.cgi.closest.filtered.bed'
    shell: '''awk '{{if (($8 > 0) && ($8 < 10)) print $1 "\t" $2 "\t" $3}}' {input} > {output}'''

rule cat_kaviar_beds:
    input: expand(DATA + 'interim/target_exp/kaviar_bed_filtered/{chr}.{{var_type}}.{{platform}}.closest.filtered.bed', chr=list(range(1,23)) + ['X','Y'])
    output: DATA + 'interim/target_exp/kaviar_bed_cat/{platform}.{var_type}.bed'
    shell: 'cat {input} > {output}'

rule cat_all_kav_beds:
    input: expand(DATA + 'interim/target_exp/kaviar_bed_cat/{platform}.{var_type}.bed', platform=('cgi'), var_type=('indel'))

rule intersect_tar_kav_beds:
    input:
        target = DATA + 'interim/target_exp/target_bed/{platform}.{var_type}.bed',
        kaviar = DATA + 'interim/target_exp/kaviar_bed_cat/{platform}.{var_type}.bed'
    output: DATA + 'interim/target_exp/target_kaviar_intersect/{platform}.{var_type}.{cov}_cov.intersect'
    shell: 'bedtools intersect -f {wildcards.cov} -r -a {input.target} -b {input.kaviar} > {output}'

rule intersect_tar_kav_all:
    input: expand(DATA + 'interim/target_exp/target_kaviar_intersect/{platform}.{var_type}.{cov}_cov.intersect', platform=('cgi'), var_type=('indel'), cov = ('0.10','0.20','0.50','0.75','0.80','0.90'))

def get_pysster_results(path):
    files_path = str(path)
    FILE_LIST = []
    for afile in os.listdir(files_path):
        FILE_LIST.append(afile)
    return FILE_LIST

rule intersect_pysster_out:
    input:
        FILE_LIST = get_pysster_results(DATA + 'interim/pysster_classification_results/'),
        pysster = DATA + 'interim/pysster_classification_results/{afile}',
        intersect = DATA + 'interim/target_exp/target_kaviar_intersect/cgi.indel.intersect' 
    output: DATA + 'interim/target_exp/target_kaviar_pysster_intersect/{afile}'
    shell: 'bedtools intersect -wa -wb -a {input.pysster} -b {input.intersect} > {output}'

rule check_intersects:
    input: expand(DATA + 'interim/target_exp/target_kaviar_pysster_intersect/{afile}', afile = FILE_LIST)
    run:
        print(FILE_LIST)

#rule intersectTargetFiles:
#    input: vcf=DATA + 'interim/target_subsets/{subset}/{subset}.vcf',
#           bed='/mnt/isilon/cbmi/variome/perry/projects/sarmadi/ahmad_exomeizer/data/interim/sql_result.status.bed'
#    output: DATA + 'interim/target_subsets/{subset}/{subset}_sql_result.status.intr'
#    shell: 'bedtools intersect -wb -a {input.vcf} -b {input.bed} > {output}'
#
#rule intersectTargetFilesAll:
#    input: expand(DATA + 'interim/target_subsets/{subset}/{subset}_sql_result.status.intr', subset=('mat_not_illumina','mat_illumina_match','mat'))
#
#rule intersectTargetDifficultByFile:
#    input: vcf=DATA + 'interim/target_subsets/{subset}/{subset}.vcf',
#           bed='/mnt/isilon/cbmi/variome/perry/projects/sarmadi/ahmad_exomeizer/data/interim/difficult/by_file/{subset2}/{subset3}.bed'
#    output: DATA + 'interim/target_subsets/{subset}/{subset}_{subset2}__{subset3}.intr'
#    shell: 'bedtools intersect -wb -a {input.vcf} -b {input.bed} > {output}'
#
#rule intersectAllTargetDifficultByFile:
#    input: expand(DATA + 'interim/target_subsets/{subset}/{subset}_{subset4}.intr', \
#                  subset=('mat_not_illumina','mat_illumina_match','mat'), \
#                  subset4=('FunctionalTechnicallyDifficultRegions__BadPromoters_gb-2013-14-5-r51-s1.sql_result.status.difficult','GCcontent__human_g1k_v37_l100_gclt25orgt65_slop50.sql_result.status.difficult','LowComplexity__AllRepeats_gt95percidentity_slop5.sql_result.status.difficult','mappability__lowmappabilityall.sql_result.status.difficult','SegmentalDuplications__segdupall.sql_result.status.difficult'))
#
#rule intersectTargetInitialDifficultFiles:
#    input:  vcf=DATA + 'interim/target_subsets/{subset}/{subset}.vcf',
#            bed='/mnt/isilon/cbmi/variome/perry/projects/sarmadi/benchmarking-tools/resources/stratification-bed-files/{subset2}/{subset3}.bed.gz'
#    output: DATA + 'interim/target_subsets/{subset}/initial_files_intersections/{subset}_{subset2}__{subset3}.intr'
#    shell:  'bedtools intersect -wb -a {input.vcf} -b {input.bed} > {output}'
#
#rule intersectAllTargetInitialDifficultFiles:
#    input: expand(DATA + 'interim/target_subsets/{subset}/initial_files_intersections/{subset}_{subset4}.intr', \
#                  subset=('mat_not_illumina','mat_illumina_match','mat'), \
#                  subset4=('mappability__lowmappabilityall','SegmentalDuplications__segdupall'))
#
#rule getTargetLineCounts:
#    input: DATA + 'interim/target_subsets/'
#    output: DATA + 'interim/file_line_counts_target'
#    shell: '''wc -l `find {input} -type f` > {output}'''
#
#rule sortTargetLineCounts:
#    input: DATA + 'interim/file_line_counts_target'
#    output: DATA + 'interim/file_line_counts_target_sorted.csv'
#    shell: 'python {SCRIPTS}line_count_sort_target.py {input} {output}'
#
