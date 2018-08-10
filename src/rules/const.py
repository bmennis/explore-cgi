import os, sys
import pandas as pd

p = os.getcwd()
if 'src' in p:
    PWD = p.split('src/')[0]
else:
    PWD = p + '/'

WORK = PWD + 'work/'
FILES = PWD + 'docs/'
SCRIPTS = PWD + 'src/scripts/'
DATA = PWD + 'data/'
LOG = PWD + 'log/'
CONFIG = PWD + 'configs/'

REGION_DIR = '/mnt/isilon/cbmi/variome/perry/projects/sarmadi/benchmarking-tools/resources/stratification-bed-files/'
GEMINI_DIR = '/mnt/isilon/cbmi/variome/bin/gemini/data/gemini_data/'

ANNO_BEDS = ['lowmappabilityall', 'notinlowmappabilityall', 'siren_similarRegions_dist1', 'segdupall', 'notinsegdupall',
             'notinrefseq_union_cds.sort', 'BadPromoters_gb-2013-14-5-r51-s1',
             'human_g1k_v37_l100_gclt30orgt55_slop50', 'human_g1k_v37_l100_gc30to55_slop50','hg19_self_chain_split.sort',
             'hg19_self_chain_split_both', 'notinAllRepeats_gt95percidentity_slop5',
             'AllRepeats_gt95percidentity_slop5', 'AllRepeats_lt51bp_gt95identity_merged',
             '1kg','20120824_combined_mask','blackTerry','dgv','dgv.short','GRCh37GenomicSuperDup.sorted','hg19.blacklist','rmsk','simpleRepeat',
             'wgEncodeCrgMapabilityAlign100mer','wgEncodeCrgMapabilityAlign24mer','wgEncodeCrgMapabilityAlign36mer',
             'wgEncodeCrgMapabilityAlign40mer','wgEncodeCrgMapabilityAlign50mer','wgEncodeCrgMapabilityAlign75mer',
             'consensusBlacklist','dukeExcludeRegions','wgEncodeDukeMapabilityUniqueness20bp','wgEncodeDukeMapabilityUniqueness35bp']
HG19_FA_NOCHR = '/mnt/isilon/cbmi/variome/reference/human/hg19/hg19NoChr.fa'
