"""Annotations for TARGET."""
rule removeTargetDuplicates:
    input: '/mnt/isilon/cbmi/variome/perry/projects/diskin/nb_convergence/data/features/mat'
    output: DATA + 'raw/mat_raw_list'
    shell:  'cat {input} | cut -f 3,4,5,6 | sort -u > {output}'

rule createTargetVcfFile:
    input: DATA + 'raw/mat_raw_list'
    output: DATA + 'raw/mat.vcf'
    shell: 'python src/scripts/vcf_create.py {input} {output}'
