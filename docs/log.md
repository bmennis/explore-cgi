###2018_06_22
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
