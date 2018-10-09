import random
import sys
import pandas as pd

kaviar_mat_file = sys.argv[1]
samp_mat_file = sys.argv[2]

dat = pd.read_csv(kaviar_mat_file, sep='\t')

random.seed(900)
samp_mat = dat.sample(8050)

samp_mat.to_csv(samp_mat_file, sep='\t')
