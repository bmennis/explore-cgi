import os, sys

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
             'notinrefseq_union_cds.sort', 'BadPromoters_gb-2013-14-5-r51-s1']
