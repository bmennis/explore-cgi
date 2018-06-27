"""Annotations for TARGET."""
rule removeTargetDuplicates:
    input: '/mnt/isilon/cbmi/variome/perry/projects/diskin/nb_convergence/data/features/mat'
    output: DATA + 'raw/mat_raw_list'
    shell:  'cat {input} | cut -f 3,4,5,6 | sort -u > {output}'

rule createTargetVcfFile:
    input: DATA + 'raw/mat_raw_list'
    output: DATA + 'raw/mat.vcf'
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
    ouput: DATA + 'raw/mat_illumina_match'
    shell: '''awk -F '\\t' '($5!="0") || ($6!="0") || ($7!="0") || ($8!="0") || ($9!="0") || ($10!="0") || ($11!="0") {{print $0}}' data/raw/mat_sorted > data/raw/mat_illumina_match'''
