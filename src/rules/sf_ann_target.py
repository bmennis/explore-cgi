"""Annotations for TARGET."""
rule removeTargetDuplicates:
    input: '/mnt/isilon/cbmi/variome/perry/projects/diskin/nb_convergence/data/features/mat'
    output: DATA + 'raw/mat_raw_list'
    shell:  'cat {input} | cut -f 3,4,5,6 | sort -u > {output}'

rule createTargetVcfFile:
    input: DATA + 'raw/mat_raw_list'
    output: DATA + 'interim/target_subsets/mat/mat.vcf'
    shell: 'python src/scripts/vcf_create.py {input} {output}'

rule filterTargetWithMatch:
    input: '/mnt/isilon/cbmi/variome/perry/projects/diskin/nb_convergence/data/features/mat'
    output: DATA + 'raw/mat_sorted'
    shell: 'cat {input} | cut -f 3,4,5,6,7,8,9,10,24,25,26 | sort -u > {output}'

rule sortTargetNotIllumina:
    input: DATA + 'raw/mat_sorted'
    output: DATA + 'raw/mat_not_illumina'
    shell: '''awk -F '\\t' '($5=="0") && ($6=="0") && ($7=="0") && ($8=="0") && ($9=="0") && ($10=="0") && ($11=="0") {{print $0}}' {input} > {output}'''

rule sortTargetIllumina:
    input: DATA + 'raw/mat_sorted'
    output: DATA + 'raw/mat_illumina_match'
    shell: '''awk -F '\\t' '($5!="0") || ($6!="0") || ($7!="0") || ($8!="0") || ($9!="0") || ($10!="0") || ($11!="0") {{print $0}}' {input} > {output}'''

rule prepTargetMatchMismatchFiles:
    input: DATA + 'raw/{subset}'
    output: DATA + 'raw/{subset}_vcf_prep'
    shell:  'cut -f 1,2,3,4 {input} > {output}'

rule prepAllTargetMatchMismatch:
    input: expand(DATA + 'raw/{subset}_vcf_prep', subset=('mat_not_illumina','mat_illumina_match'))

rule createTargetMatchMismatchVcfFile:
    input: DATA + 'raw/{subset}_vcf_prep'
    output: DATA + 'interim/target_subsets/{subset}/{subset}.vcf'
    shell: 'python src/scripts/vcf_create.py {input} {output}'

rule createAllMatchMismatchVcfFiles:
    input: expand(DATA + 'interim/target_subsets/{subset}/{subset}.vcf', subset=('mat_not_illumina','mat_illumina_match')) 

rule intersectFiles:
    input: vcf=DATA + 'interim/target_subsets/{subset}/{subset}.vcf',
           bed='/mnt/isilon/cbmi/variome/perry/projects/sarmadi/ahmad_exomeizer/data/interim/sql_result.status.bed'
    output: DATA + 'interim/target_subsets/{subset}/{subset}_hgmd.intr'
    shell: 'bedtools intersect -wb -a {input.vcf} -b {input.bed} > {output}'

rule intersectFilesAll:
    input: expand(DATA + 'interim/target_subsets/{subset}/{subset}_hgmd.intr', subset=('mat_not_illumina','mat_illumina_match','mat'))

rule intersectDifficultByFile:
    input: vcf=DATA + 'interim/target_subsets/{subset}/{subset}.vcf',
         bed='/mnt/isilon/cbmi/variome/perry/projects/sarmadi/ahmad_exomeizer/data/interim/difficult/by_file/{subset2}/{subset3}.bed'
    output: DATA + 'interim/target_subsets/{subset}/{subset2}_{subset3}.intr'
    shell: 'bedtools intersect -wb -a {input.vcf} -b {input.bed} > {output}'

rule intersectAllDifficultByFile:
    input: expand(DATA + 'interim/target_subsets/{subset}/{subset2}_{subset3}.intr', \
                  subset=('mat_not_illumina','mat_illumina_match','mat'),
                  subset2=('FunctionalTechnicallyDifficultRegions','GCcontent','LowComplexity','mappability','SegmentalDuplications'),
                  subset3=('BadPromoters_gb-2013-14-5-r51-s1.sql_result.status.difficult','human_g1k_v37_l100_gclt25orgt65_slop50.sql_result.status.difficult','AllRepeats_gt95percidentity_slop5.sql_result.status.difficult','lowmappabilityall.sql_result.status.difficult','segdupall.sql_result.status.difficult'))

rule intersectInitialDifficultFiles:
    input:  vcf=DATA + 'interim/target_subsets/{subset}/{subset}.vcf',
            bed='/mnt/isilon/cbmi/variome/perry/projects/sarmadi/benchmarking-tools/resources/stratification-bed-files/{subset2}/{subset3}.bed.gz'
    output: DATA + 'interim/target_subsets/{subset}/initial_files_intersections/{subset}_{subset3}.intr'
    shell:  'bedtools intersect -wb -a {input.vcf} -b {input.bed} > {output}'

rule all:
    input: expand(DATA + 'interim/target_subsets/{subset}/initial_files_intersections/{subset}_{subset3}.intr', \
                  subset=('mat_not_illumina','mat_illumina_match','mat'), 
                  subset2=('mappability','SegmentalDuplications'), 
                  subset3=('lowmappabilityall','segdupall'))
