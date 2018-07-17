"""Annotations for kaviar."""

rule unzip_kaviar:
    input:  '/mnt/isilon/cbmi/variome/bin/gemini/data/gemini_data/Kaviar-160204-Public-hg19.vt.vcf.gz'
    output: DATA + 'interim/kaviar.vcf'
    shell:  'gunzip -c {input} > {output}'

rule split_kaviar_chrom:
    input:  i = DATA + 'interim/kaviar.vcf'
    output: expand(DATA + 'interim/kaviar/{chrom}.full.vcf', chrom=list(range(1,23)) + ['X', 'Y', 'M'])
    run:
        outs = {}
        for chrom in list(range(1,23)) + ['X', 'Y', 'M']:
            outs[str(chrom)] = open(DATA + 'interim/kaviar/' + str(chrom) + '.full.vcf', 'w')
        with open(input.i) as f:
            for line in f:
                if line[0] == '#':
                    for chrom in outs:
                        print(line.strip(), file=outs[chrom])
                else:
                    chrom = line.split('\t')[0]
                    print(line.strip(), file=outs[chrom])

rule mk_short_kaviar:
    input:  i = DATA + 'interim/kaviar/{chr}.full.vcf'
    output: o = DATA + 'interim/kaviar/{chr}.short.vcf'
    run:
        with open(input.i) as f, open(output.o, 'w') as fout:
            i = 0
            for line in f:
                print(line.strip(), file=fout)
                i += 1
                if i == 5000000:
                    break

rule kaviar_vcfanno:
    input:   vcf = DATA + 'interim/kaviar/{chr}.{set}.vcf',
             config = CONFIG + 'kaviar_vcfanno.conf'
    output:  DATA + 'interim/kaviar_anno/{chr}.{set}.vcf'
    threads: 3
    shell:   'vcfanno -p {threads} -base-path {GEMINI_DIR} {input.config} {input.vcf} > {output}'

rule parse_kaviar_anno:
    input:  DATA + 'interim/kaviar_anno/{chr}.{set}.vcf'
    output: DATA + 'interim/kaviar_anno_parsed/{chr}.{set}.mat'
    shell:  'python {SCRIPTS}mk_kaviar_matrix.py {input} {output}'

rule mk_kaviar_bed:
    input:  i = DATA + 'interim/kaviar_anno_parsed/{chr}.{set}.mat'
    output: cgi = DATA + 'interim/kaviar_bed/{chr}.{set}.cgi.bed',
            both = DATA + 'interim/kaviar_bed/{chr}.{set}.both.bed',
    run:
        dat = pd.read_csv(input.i, sep='\t')
        features = [x for x in dat.columns if not x in ('kaviar_status', 'var_type', 'chrom', 'pos', 'ref', 'alt')]
        dat = dat.drop_duplicates(subset=['chrom', 'pos'])
        dat_ahmad = dat[((dat['ahmad_status'] == 1) | (dat['ahmad_status'] == 0)) ]
        dat_ahmad.loc[:, 'pos_minus'] = dat_ahmad['pos'] - 150
        dat_ahmad.loc[:, 'pos_plus'] = dat_ahmad['pos'] + 150
        cgi_len = len(dat_ahmad[dat_ahmad.kaviar_status=='cgi'])
        both_len = len(dat_ahmad[dat_ahmad.kaviar_status=='both'])
        size = min((cgi_len, both_len))
        dat_ahmad[dat_ahmad.kaviar_status=='cgi'].sample(size)[['chrom', 'pos_minus', 'pos_plus']].to_csv(output.cgi, index=False, header=False, sep='\t')
        dat_ahmad[dat_ahmad.kaviar_status=='both'].sample(size)[['chrom', 'pos_minus', 'pos_plus']].to_csv(output.both, index=False, header=False, sep='\t')

rule mk_kaviar_fasta:
    input:  bed = DATA + 'interim/kaviar_bed/{chr}.{set}.bed',
            fa = HG19_FA_NOCHR
    output: DATA + 'interim/kaviar_fa/{chr}.{set}.fa'
    shell:  'bedtools getfasta -fi {input.fa} -bed {input.bed} > {output}'

rule upper_fa:
    input:  i = DATA + 'interim/kaviar_fa/{aset}.fa'
    output: o = DATA + 'interim/kaviar_fa_upper/{aset}.fa'
    run:
        with open(input.i) as f, open(output.o, 'w') as fout:
            for line in f:
                if line[0] == '>':
                    print(line.strip(), file=fout)
                else:
                    print(line.upper().strip(), file=fout)

rule collapse_kaviar_fasta:
    input:  expand(DATA + 'interim/kaviar_fa_upper/{chr}.{{aset}}.fa', chr=list(range(1,23)) + ['X', 'Y', 'M'])
    output: DATA + 'interim/kaviar_fa_gz/{aset}.fa.gz'
    shell:  'cat {input} | gzip - > {output}'

rule all_kaviar_fa:
    input: expand(DATA + 'interim/kaviar_fa_gz/{aset}.fa.gz', aset=('short.both', 'short.cgi'))

rule mk_short_kaviars:
    input:  expand(DATA + 'interim/kaviar_anno_parsed/{chr}.short.mat', chr=list(range(1,23)) + ['X', 'Y', 'M'])
    output: o = DATA + 'interim/kaviar.mat'
    run:
        with open(list(input)[0]) as f:
            header = f.readline().strip()
        with open(output.o, 'w') as fout:
            print(header, file=fout)
            for i in input:
                with open(i) as f:
                    f.readline()
                    for line in f:
                        print(line.strip(), file=fout)

rule split_kaviar:
    input:  DATA + 'interim/kaviar.vcf'
    output: expand(DATA + 'interim/kaviar_subsets/{subset}.vcf', subset=('cgi_only_sources', 'illumina_only_sources', 'cgi_illumina_sources','no_sources'))
    shell:  'python {SCRIPTS}vcf_sort.py {input} {DATA}interim/kaviar_subsets/'

rule filterIlluminaOnlySources:
    input:  DATA + 'interim/kaviar_subsets/illumina_only_sources.vcf'
    output: DATA + 'interim/kaviar_subsets/illumina_only_sources_filtered.vcf'
    shell:  '''awk -F'\\t' -vOFS='\\t' '{{ if ($5 ~ /^<.*/) $5="."}}1' {input} > {output}'''

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

