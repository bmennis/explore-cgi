"""Annotations for kaviar."""

rule unzip_kaviar:
    input:  '/mnt/isilon/cbmi/variome/bin/gemini/data/gemini_data/Kaviar-160204-Public-hg19.vt.vcf.gz'
    output: DATA + 'interim/kaviar.vcf'
    shell:  'gunzip -c {input} > {output}'

rule split_kaviar:
    input:  DATA + 'interim/kaviar.vcf'
    output: expand(DATA + 'interim/kaviar_subsets/{subset}.vcf', subset=('cgi_only_sources', 'illumina_only_sources', 'cgi_illumina_sources','no_sources'))
    shell:  'python {SCRIPTS}vcf_sort.py {input} {DATA}interim/kaviar_subsets/'

rule filterIlluminaOnlySources:
    input: DATA + 'interim/kaviar_subsets/illumina_only_sources.vcf'
    output: DATA + 'interim/kaviar_subsets/illumina_only_sources_filtered.vcf'
    shell: '''awk -F'\\t' -vOFS='\\t' '{{ if ($5 ~ /^<.*/) $5="."}}1' {input} > {output}'''

rule intersectInitialDifficultFiles:
    input:  vcf=DATA + 'interim/kaviar_subsets/{subset}.vcf',
            bed='/mnt/isilon/cbmi/variome/perry/projects/sarmadi/benchmarking-tools/resources/stratification-bed-files/{subset2}/{subset3}.bed.gz'
    output: DATA + 'interim/kaviar_subsets/{subset}/initial_files_intersections/{subset}_{subset2}__{subset3}.intr'
    shell:  'bedtools intersect -wb -a {input.vcf} -b {input.bed} > {output}'

rule all:
    input: expand(DATA + 'interim/kaviar_subsets/{subset}/initial_files_intersections/{subset}_{subset4}.intr', \
                  subset=('cgi_only_sources','illumina_only_sources_filtered','cgi_illumina_sources','no_sources'), \
                  subset4=('mappability__lowmappabilityall','SegmentalDuplications__segdupall'))

rule intersectFiles:
    input: vcf=DATA + 'interim/kaviar_subsets/{subset}.vcf',
           bed='/mnt/isilon/cbmi/variome/perry/projects/sarmadi/ahmad_exomeizer/data/interim/sql_result.status.bed'
    output: DATA + 'interim/kaviar_subsets/{subset}/{subset}_sql_result.status.intr'
    shell: 'bedtools intersect -wb -a {input.vcf} -b {input.bed} > {output}'

rule intersectFilesAll:
    input: expand(DATA + 'interim/kaviar_subsets/{subset}/{subset}_sql_result.status.intr', \
                  subset=('cgi_only_sources','illumina_only_sources_filtered','cgi_illumina_sources','no_sources'))

rule intersectDifficultByFile:
    input: vcf=DATA + 'interim/kaviar_subsets/{subset}.vcf',
           bed='/mnt/isilon/cbmi/variome/perry/projects/sarmadi/ahmad_exomeizer/data/interim/difficult/by_file/{subset2}/{subset3}.bed'
    output: DATA + 'interim/kaviar_subsets/{subset}/{subset}_{subset2}__{subset3}.intr'
    shell: 'bedtools intersect -wb -a {input.vcf} -b {input.bed} > {output}'

rule intersectAllDifficultByFile:
    input: expand(DATA + 'interim/kaviar_subsets/{subset}/{subset}_{subset4}.intr', \
                  subset=('cgi_only_sources','illumina_only_sources_filtered','cgi_illumina_sources','no_sources'), \
                  subset4=('FunctionalTechnicallyDifficultRegions__BadPromoters_gb-2013-14-5-r51-s1.sql_result.status.difficult','GCcontent__human_g1k_v37_l100_gclt25orgt65_slop50.sql_result.status.difficult','LowComplexity__AllRepeats_gt95percidentity_slop5.sql_result.status.difficult','mappability__lowmappabilityall.sql_result.status.difficult','SegmentalDuplications__segdupall.sql_result.status.difficult'))

rule sortAlleleFrequency:
    input:  expand(DATA + 'interim/kaviar_subsets/{subset}.vcf', subset=('cgi_only_sources','cgi_illumina_sources','illumina_only_sources_filtered','no_sources'))
    output: DATA + 'interim/Kaviar_subset_allele_frq_complete.csv'
    run:
        for afile in input:
            shell('python {SCRIPTS}vcf_sort_allele.py {afile} {output}')

rule getKaviarVcfLineCounts:
    input: DATA + 'interim/kaviar.vcf'
    output: DATA + 'interim/kaviar_vcf_line_counts'
    shell: 'wc -l {input} > {output}'

rule getKaviarFileLineCounts:
    input: DATA + 'interim/kaviar_subsets/'
    output: DATA + 'interim/kaviar_file_line_counts'
    shell: '''wc -l `find {input} -type f` > {output}'''

rule getKaviarTotalLineCounts:
    input: txt1=DATA + 'interim/kaviar_vcf_line_counts',
           txt2=DATA + 'interim/kaviar_file_line_counts'
    output: DATA + 'interim/kaviar_total_line_counts'
    shell: 'cat {input.txt1} {input.txt2} > {output}'

rule sortKaviarLineCounts:
    input: DATA + 'interim/kaviar_total_line_counts'
    output: DATA + 'interim/kaviar_line_counts_sorted.csv'
    shell: 'python {SCRIPTS}line_count_sort.py {input} {output}'

# rule sortAlleleFrequencyKaviar:
#    input: vcf=DATA + '/mnt/isilon/cbmi/variome/bin/gemini/data/gemini_data/Kaviar-160204-Public-hg19.vt.vcf.gz',
#    output: DATA + 'raw/Kaviar_allele_frq.csv'
#    shell: 'python src/scripts/vcf_sort_allele.py {input.vcf} {output.csv}'

