"""Annotations for kaviar."""

rule unzip_kaviar:
    input:  '/mnt/isilon/cbmi/variome/bin/gemini/data/gemini_data/Kaviar-160204-Public-hg19.vt.vcf.gz'
    output: DATA + 'interim/kaviar.vcf'
    shell:  'gunzip -c {input} > {output}'

rule split_kaviar:
    input:  DATA + 'interim/kaviar.vcf'
    output: expand(DATA + 'interim/kaviar_subsets/{subset}.vcf', subset=('cgi_only_sources', 'illumina_only_sources', 'cgi_illumina_sources'))
    shell:  'python {SCRIPTS}vcf_sort.py {input} {DATA}interim/kaviar_subsets/'

rule filterIlluminaOnlySources:
    input: DATA + 'interim/kaviar_subsets/illumina_only_sources.vcf'
    output: DATA + 'interim/kaviar_subsets/illumina_only_sources.vcf'
    shell: '''awk -F'\\t' -vOFS='\\t' '{ if ($5 ~ /^<.*/) $5="."}1' {input} > tmp && mv tmp {output}'''

rule intersectInitialDifficultFiles:
    input:  vcf=DATA + 'interim/kaviar_subsets/{subset}.vcf',
            bed='/mnt/isilon/cbmi/variome/perry/projects/sarmadi/benchmarking-tools/resources/stratification-bed-files/{subset2}/{subset3}.bed.gz'
    output: DATA + 'interim/kaviar_subsets/{subset}/initial_files_intersections/{subset}_{subset3}.intr'
    shell:  'bedtools intersect -wb -a {input.vcf} -b {input.bed} > {output}'

# rule all:
#     input: expand(DATA + 'interim/{subset}/initial_files_intersections/{subset}_{subset3}.intr', subset=('cgi_only_sources','illumina_only_sources','cgi_illumina_sources','illumina_only_sources_test','mat'), subset2=('mappability','SegmentalDuplications'), subset3=('lowmappabilityall','segdupall'))

rule intersectFiles:
    input: vcf=DATA + 'interim/kaviar_subsets/{subset}.vcf',
           bed='/mnt/isilon/cbmi/variome/perry/projects/sarmadi/ahmad_exomeizer/data/interim/sql_result.status.bed'
    output: DATA + 'interim/kaviar_subsets/{subset}/{subset}_hgmd.intr'
    shell: 'bedtools intersect -wb -a {input.vcf} -b {input.bed} > {output}'

rule intersectFilesAll:
    input: expand(DATA + 'interim/{subset}/{subset}_hgmd.intr', subset=('cgi_only_sources','illumina_only_sources_filtered','cgi_illumina_sources'))

# rule intersectDifficultByFile:
#     input: vcf=DATA + 'raw/{subset}.vcf',
# 	   bed='/mnt/isilon/cbmi/variome/perry/projects/sarmadi/ahmad_exomeizer/data/interim/difficult/by_file/{subset2}/{subset3}.bed'
#     output: DATA + 'interim/{subset}/{subset2}_{subset3}.intr'
#     shell: 'bedtools intersect -wb -a {input.vcf} -b {input.bed} > {output}'

# rule intersectAllDifficultByFile:
#     input: expand(DATA + 'interim/{subset}/{subset}_{subset3}.intr', \
#                   subset=('cgi_only_sources','illumina_only_sources','illumina_only_sources_test','cgi_illumina_sources','mat'),
#                   subset2=('FunctionalTechnicallyDifficultRegions','GCcontent','LowComplexity','mappability','SegmentalDuplications'),
#                   subset3=('BadPromoters_gb-2013-14-5-r51-s1.sql_result.status.difficult','human_g1k_v37_l100_gclt25orgt65_slop50.sql_result.status.difficult','AllRepeats_gt95percidentity_slop5.sql_result.status.difficult','lowmappabilityall.sql_result.status.difficult','segdupall.sql_result.status.difficult'))

rule sortAlleleFrequency:
    input:  expand(DATA + 'interim/kaviar_subsets/{subset}.vcf', subset=('cgi_only_sources','cgi_illumina_sources','illumina_only_sources'))
    output: DATA + 'interim/Kaviar_subset_allele_frq.csv'
    run:
        for afile in input:
            shell('python {SCRIPTS}vcf_sort_allele.py {afile} {output}')

# rule sortAlleleFrequencyAll:
#     input: expand(DATA + 'raw/Kaviar_subset_allele_frq.csv', subset=('cgi_only_sources','cgi_illumina_sources','illumina_only_sources')

# rule sortAlleleFrequencyKaviar:
#     input: vcf=DATA + '/mnt/isilon/cbmi/variome/bin/gemini/data/gemini_data/Kaviar-160204-Public-hg19.vt.vcf.gz',
#     output: DATA + 'raw/Kaviar_allele_frq.csv'
#     shell: 'python src/scripts/vcf_sort_allele.py {input.vcf} {output.csv}'

