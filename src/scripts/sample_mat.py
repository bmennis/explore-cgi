import random
import sys
import pandas as pd

kaviar_mat_file = sys.argv[1]
samp_mat_file = sys.argv[2]
unsampled_mat_file = sys.argv[3]

cols = ['chrom', 'pos', 'ref', 'alt']
dat = pd.read_csv(kaviar_mat_file, sep='\t')
random.seed(900)
samp_mat = dat.sample(8050)

m = pd.merge(dat, samp_mat[cols], on=cols, how='left', indicator=True)
m[m._merge=='left_only'].drop(['_merge'], axis=1).to_csv(unsampled_mat_file, index=False, sep='\t')
samp_mat.to_csv(samp_mat_file, sep='\t')
