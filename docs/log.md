### 2018_11_8
* checked allele frequencies of indels after filtering and vast majority very low frq
* looking into Target dataset
* after filtering out duplicates get counts for indels and snvs
* go back to Kaviar set and get allele frequencies for snvs of cgi only data, want to filter snvs for good ahmad status
* With both counts look into possibly seeing if we can get a merge of the two data sets and find common reads

### 2018_10_24
* we recognized that training of the pysster model would need downsampling of the unverified cgi data set to roughly the verified set
* once the data was downsampled, the pysster model was trained on all additional features to this point and saved
* ran the pysster model on the unsampled data for validation.
* to explore the performance on differing features and annotations we ran pysster with no added features, all added features and only reference sequence
* pysster appears to remove 'useless' (all flagged 0 or 1) features when training a model, beware of errors running saved model for different amounts of features
* need to investigate the frequency cutoff limit

### 2018_08_28
* ran classifier filtering cgi indels out that are between 0 and 10bp from both or ill indels
* need to run classifier on full set of indel variants with no filtering
* need to look into g quadruplex annotation for added feature

### 2018_07_27
* need to split indels and snvs
* why are snvs so hard? see if they fix the indel pattern

### 2018_07_26
* run w/ additional features
```
qsub -q gpu.q -l GPU=1 -l h_vmem=50G -l mem_free=50G -e ./perry_err -o ./perry_out pysster_script.sh
```

### 2018_07_19
* pysster is slowing down
* trying pysster run to max pr during grid search and focus on only positions w/ good ahmad coverage

### 2018_07_18
* tried pysster w/ nuc alphabet, pyrimidine vs purine, and cg vs at, and 300 nuc window, and 100 kernels. got 73% accuracy. Ran 14 hrs.
* Can request 4 gpus. Trying to maximize precision of cgi only for next run. Need repeat features.

### 2018_07_17
* tried pysster w/ nuc alphabet, 150 nuc window, and 20 as kernel length. got 73% accuracy. Ran 3 hrs.

### 2018_07_13
* no good feature combo splits cgi only from cgi and illumina
* ahmad's good/poor regions would work well, but there are not enough cgi only variants in the poor regions
* going to try deep learning: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty222/4962494
* limiting af helps. af > 2^-13 gives 1800 negs, and 73% accuracy. ahmad's region is necessary for the boost, so we are seeing the results of a lack of illumina coverage.
* af > 2^-14 gives 3500 negs (inside ahmad's regions) and at 76% accuracy
* so cgi only at low and high af behave diferently
* need to check low freq variants. how many are there af < 2**-14? 30K cgi only and 8.5K both. Accuracy is now 80%.
* but putting low and high together results in an accuracy of .5? No that accuracy is .88, so something was off in the prior runs.

### 2018_07_12
* short runs for kaviar matrix. I see no subs, which worries me.
* probably want to limit by position. Yes many variants are tied to the same position.

### 2018_07_11
* Kaviar has updated the sources so variants that had no source now have illumina or cgi sources. Will ignore for now. Ex chr1:877827
* [DL examples](https://github.com/hussius/deeplearning-biology#genomics_variant-calling)

### 2018_07_09
* Finish script to sort line counts of kaviar subsets
* Modify kaviar line count script to apply to target data
* Create rules in respective snake files to get line counts from kaviar and target files created in pipeline
* Create rules to sort line counts for kaviar and target data and put into csv files
* Pull lines from kaviar.vcf without data source which appear to be the lines missing from subset files.

### 2018_06_29
* Finish running intersection rules from snakemake files for kaviar files
* Run sortAlleleFrequency rule on Kaviar data to generate allele frequency file
* Next need to plot allele frequencies when finished
* Also need counts of variants in subsets

### 2018_06_28
*  Worked on snakemake files for kaviar and target annotation
*  Updated intersection rules to properly recognize directories and files and add those names to the outputs of the intersections
*  Tested the snakemake files and rules and fix errors
*  Run rules from snakemake files to process correct subsets of kaviar files

### 2018_06_26
* Added Target filtering and sorting rules to sf_ann_target.py to filter illumina matches and mismatches as well as sort them to different files.
* Work on updating Kaviar annotation pipeline to include filtering illumina data source variant notation.  This notation causes issues with bedtools intersections and working on recreating the command to accurately and efficiently find and replace that notation.

### 2018_06_22
* Used vcf sort allele script to pull source information for allele frequencies and allele frequencies of each read from each source vcf file.
* Used allele plot R script to plot allele sources vs frequencies for cgi only, cgi and illumina, and illumina only.
* Updated Snakemake file to include allele sorting, allele plotting, along with variant source sorting.

### 2018_06_18
* Fixed error in vcf file (blank line above header) run bedtools to annotate Target data
* Intersect Target data with bed files of good and poor mapping regions and bed files without those regions mapped

### 2018_06_15
* Continue to work on Target vcf file
* Work on graphing allele frequencies with ggplot

### 2018_06_14
* Create vcf file of Target data
* bedtools not running on vcf file, need to figure out what is causing errors

### 2018_06_13
* Create Snakemake file for project pipeline
* Add rules to Snakemake file for bedfile intersections of platform source vcf files to bed files via bedtools
* Parse Target data and create vcf file
* Need to annotate Target vcf like Kaviar platform sources and compare Target data to Kaviar
* Need to graph allele frequencies of Kaviar data and subsets and compare

### 2018_06_12
* Intersect enumerated sources of Kaviar files with multiple bed files using bedtools to annotate source variants for poor/good regions
* Removed ALT from illumina only sources of notation <CNx> and replaced with . due to errors from intersect runs
* Intersect enumerated sources of Kaviar files with multiple difficult region bed files for annotation without poor/good regions 
* Counted the number of variants of each enumerated source type both in and out of region of each bedfile without poor/good regions mapped

### 2018_06_11
* Review script and conditions of enumeration
* Run script on vcf file again on cluster

### 2018_06_08
* Run script for enumeration of kaviar sources on vcf file

### 2018_05_05
* First day - setup repo and conda env
* Start w/ enumeratoin of kaviar sources
