"""Organize region files for vcfanno.
/mnt/isilon/cbmi/variome/perry/projects/sarmadi/benchmarking-tools/resources/stratification-bed-files
"""

rule prep_low_mappability:
    input:  REGION_DIR + 'mappability/{bed}.bed.gz'
    output: GEMINI_DIR + '{bed,lowmappabilityall|notinlowmappabilityall|siren_similarRegions_dist1}.bed.gz'
    shell:  'cp {input} {output}'

rule prep_functional_regions:
    input:  REGION_DIR + 'FunctionalRegions/{bed}.bed.gz'
    output: GEMINI_DIR + '{bed,refseq_union_cds.sort|notinrefseq_union_cds.sort}.bed.gz'
    shell:  'cp {input} {output}'

rule prep_functional_technically_difficult_regions:
    input:  REGION_DIR + 'FunctionalTechnicallyDifficultRegions/{bed}.bed.gz'
    output: GEMINI_DIR + '{bed,BadPromoters_gb-2013-14-5-r51-s1}.bed.gz'
    shell:  'cp {input} {output}'

rule prep_segmental_duplications:
    input: REGION_DIR + 'SegmentalDuplications/{bed}.bed.gz'
    output: GEMINI_DIR + '{bed,hg19_self_chain_split.sort|segdupall|notinsegdupall|hg19_self_chain_split_both}.bed.gz'
    shell: 'cp {input} {output}'

rule prep_gc_content:
    input: REGION_DIR + 'GCcontent/{bed}.bed.gz'
    output: GEMINI_DIR + '{bed,human_g1k_v37_l100_gclt30orgt55_slop50|human_g1k_v37_l100_gc30to55_slop50}.bed.gz'
    shell: 'cp {input} {output}'

rule prep_low_complexity:
    input: REGION_DIR + 'LowComplexity/{bed}.bed.gz'
    output: GEMINI_DIR + '{bed,AllRepeats_gt95percidentity_slop5|notinAllRepeats_gt95percidentity_slop5|AllRepeats_lt51bp_gt95identity_merged|SimpleRepeat_homopolymer_6to10|SimpleRepeat_homopolymer_gt10|SimpleRepeat_imperfecthomopolgt10_slop5}.bed.gz'
    shell: 'cp {input} {output}'

rule prep_exomizer:
    input:  '/home/evansj/me/projects/sarmadi/ahmad_exomeizer/data/interim/sql_result.status.bed'
    output: GEMINI_DIR + 'ahmad_exomizer.bed.gz'
    shell:  'cut -f 1-3,18,19 {input} | bgzip > {output}'

rule prep_blacklist_genome:
    input: '/home/ennisb/me/blacklist_genome/{bed}.bed'
    output: GEMINI_DIR + '{bed,1kg|20120824_combined_mask|blackTerry|dgv|dgv.short|GRCh37GenomicSuperDup.sorted|hg19.blacklist|rmsk|simpleRepeat}.bed.gz'
    shell: 'bgzip -c {input} > {output}'

rule prep_alignability:
    input: DATA + 'interim/alignability/{bed}.bed'
    output: GEMINI_DIR + '{bed,wgEncodeCrgMapabilityAlign100mer|wgEncodeCrgMapabilityAlign24mer|wgEncodeCrgMapabilityAlign36mer|wgEncodeCrgMapabilityAlign40mer|wgEncodeCrgMapabilityAlign50mer|wgEncodeCrgMapabilityAlign75mer}.bed.gz'
    shell: 'bgzip -c {input} > {output}'

rule prep_excludable:
    input: DATA + 'interim/excludable/{bed}.bed'
    output: GEMINI_DIR + '{bed,consensusBlacklist|dukeExcludeRegions}.bed.gz'
    shell: 'bgzip -c {input} > {output}'

rule prep_uniqueness:
    input: DATA + 'interim/uniqueness/{bed}.bed'
    output: GEMINI_DIR + '{bed,wgEncodeDukeMapabilityUniqueness20bp|wgEncodeDukeMapabilityUniqueness35bp}.bed.gz'
    shell: 'bgzip -c {input} > {output}'

rule prep_gquad:
    input: DATA + 'interim/gquad_bed/{bed}.bed'
    output: GEMINI_DIR + '{bed,Na_PDS_plus_hits_intersect|Na_K_plus_hits_intersect|GSE63874_Na_PDS_minus_hits_intersect|GSE63874_Na_K_minus_hits_intersect}.bed.gz'
    shell: 'bgzip -c {input} > {output}'

rule prep_add_repeats:
    input: DATA + 'interim/add_repeat_feat/{bed}.bed'
    output: GEMINI_DIR + '{bed,interrupted_repeats|microsat|cpg_island_ext}.bed.gz'
    shell: 'bgzip -c {input} > {output}'

rule prep_microsat_list:
    input: GEMINI_DIR + 'hg19.2014.regions.bed'
    output: GEMINI_DIR + 'hg19.2014.regions.bed.gz'
    shell: 'bgzip -c {input} > {output}'

rule tabix_regions:
    """Index any bed file"""
    input:  GEMINI_DIR + '{bedfile}.bed.gz'
    output: GEMINI_DIR + '{bedfile}.bed.gz.tbi'
    shell:  'tabix -p bed {input}'

rule all_vcfanno_files:
    input: GEMINI_DIR + 'ahmad_exomizer.bed.gz.tbi',
           files=expand(GEMINI_DIR + '{bedfile}.bed.gz.tbi', bedfile=ANNO_BEDS),
    output: o = CONFIG + 'kaviar_vcfanno.conf'
    run:
        with open(output.o, 'w') as fout:
            anno_entry = """[[annotation]]
file="ahmad_exomizer.bed.gz"
names=["ahmad_region","ahmad_region_mq"]
columns=[4,5]
ops=["self","self"]
"""

            microsat_entry = """[[annotation]]
file="hg19.2014.regions.bed.gz"
names=["microsat_info"]
columns=[4]
ops=["self"]
"""
            print(anno_entry, file=fout)
            print(microsat_entry, file=fout)
            for afile in input.files:
                file_name = afile.split('/')[-1].strip('.tbi')
                name = file_name.split('.bed')[0]
                anno_entry = """[[annotation]]
file="%s"
names=["%s"]
columns=[2]
ops=["flag"]
""" % (file_name, name)
                print(anno_entry, file=fout)
