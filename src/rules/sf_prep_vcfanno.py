"""Organize region files for vcfanno.
/mnt/isilon/cbmi/variome/perry/projects/sarmadi/benchmarking-tools/resources/stratification-bed-files
"""

rule low_mappability:
    input: REGION_DIR + 'mappability/{bed}.bed.gz'
    output: GEMINI_DIR + '{bed,lowmappabilityall|notinlowmappabilityall|siren_similarRegions_dist1}.bed.gz'
    shell: 'cp {input} {output}'

rule tabix_regions:
    """Index any bed file"""
    input:  GEMINI_DIR + '{bedfile}.bed.gz'
    output: GEMINI_DIR + '{bedfile}.bed.gz.tbi'
    shell:  'tabix -p bed {input}'

ANNO_BEDS = ('lowmappabilityall', 'notinlowmappabilityall', 'siren_similarRegions_dist1')
rule all_vcfanno_files:
    input:  expand(GEMINI_DIR + '{bedfile}.bed.gz.tbi', bedfile=ANNO_BEDS)
    output: o = CONFIG + 'kaviar_vcfanno.conf'
    run:
        with open(output.o, 'w') as fout:
            for afile in input:
                file_name = afile.split('/')[-1].strip('.tbi')
                name = file_name.split('.bed')[0]
                anno_entry = """[[annotation]]
file="%s"
names=["%s"]
columns=[2]
ops=["flag"]
""" % (file_name, name)
                print(anno_entry, file=fout)
