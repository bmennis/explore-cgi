"""Organize region files for vcfanno.
/mnt/isilon/cbmi/variome/perry/projects/sarmadi/benchmarking-tools/resources/stratification-bed-files
"""

rule prep_low_mappability:
    input:  REGION_DIR + 'mappability/{bed}.bed.gz'
    output: GEMINI_DIR + '{bed,lowmappabilityall|notinlowmappabilityall|siren_similarRegions_dist1}.bed.gz'
    shell:  'cp {input} {output}'

rule prep_functional_regions:
    input: REGION_DIR + 'FunctionalRegions/{bed}.bed.gz'
    output: GEMINI_DIR + '{bed,refseq_union_cds.sort|notinrefseq_union_cds.sort}.bed.gz'
    shell: 'cp {input} {output}'

rule prep_functional_technically_difficult_regions:
    input: REGION_DIR + 'FunctionalTechnicallyDifficultRegions/{bed}.bed.gz'
    output: GEMINI_DIR + '{bed,BadPromoters_gb-2013-14-5-r51-s1}.bed.gz'
    shell: 'cp {input} {output}'

rule prep_segmental_duplications:
    input: REGION_DIR + 'SegmentalDuplications/{bed}.bed.gz'
    ouput: GEMINI_DIR + '{bed,hg19_self_chain_split.sort|segdupall|notinsegdupall|hg19_self_chain_split_both}.bed.gz'
    shell: 'cp {input} {output}'

#Working on this rule
#rule prep_gc_content:
#    input: REGION_DIR + 'GCcontent/{bed}.bed.gz'
#    output: GEMINI_DIR + '{bed,

rule tabix_regions:
    """Index any bed file"""
    input:  GEMINI_DIR + '{bedfile}.bed.gz'
    output: GEMINI_DIR + '{bedfile}.bed.gz.tbi'
    shell:  'tabix -p bed {input}'

ANNO_BEDS = ('lowmappabilityall', 'notinlowmappabilityall', 'siren_similarRegions_dist1', 'refseq_union_cds.sort', 'notinrefseq_union_cds.sort', 'BadPromoters_gb-2013-14-5-r51-s1')
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
