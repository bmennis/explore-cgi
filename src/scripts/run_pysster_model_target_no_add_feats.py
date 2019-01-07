import os
from pysster.Data import Data
from pysster import utils
from IPython.display import Image
import pandas as pd
import gzip
DATA = "/mnt/isilon/dbhi_bfx/perry/brian/"

def get_fa_info(fa_files):
    ###Function to read in fa file and take all of the chrom, pos start, pos end info for classification results bed file###
    chroms = []
    starts = []
    ends = []
    for a_file in fa_files:
        with gzip.open(a_file, 'rt') as inf:
            for line in inf:
                line = line.strip()
                if line[0] == '>':
                    chrom = line.split('>')[1][0]
                    pos_start = line.split(':')[1].split('-')[0]
                    pos_end = line.split(':')[1].split('-')[1]
                    chroms.append(chrom)
                    starts.append(pos_start)
                    ends.append(pos_end)
                else:
                    pass

    return chroms, starts, ends

#establish output directory and take output directory name for csv file of labels and predictions
output_folder =  DATA + "explore-cgi/data/interim/cgi_ind_exp/pysster_output/model_run_12_21_18_kav_cgi_8k_samp_kav_both_8k_samp/"
if not os.path.isdir(output_folder):
    os.makedirs(output_folder)

class_file_name = output_folder.split('/')[-2]

#load the pysster prediction model
model = utils.load_model("/mnt/isilon/dbhi_bfx/perry/brian/explore_cgi/data/interim/target_exp/pysster_output/target_train_model_no_add_feats_12_21_18/model.pkl")


#load the dataset as data
data = Data([DATA + "explore-cgi/data/interim/cgi_ind_exp/pysster_fa/cgi.indel.sample.fa.gz", DATA + "explore-cgi/data/interim/cgi_ind_exp/pysster_fa/both.indel.sample.fa.gz"], ("ACGT", "XDI"))


#run the model of pysster on all of the data set
predictions = model.predict(data, "all")

labels = data.get_labels("all")

#create a data frame of the labels and predictions for each variant and write that dataframe to a csv file for analysis
#fa_files = [DATA + "explore-cgi/data/interim/cgi_ind_exp/pysster_fa/cgi.indel.sample.fa.gz", DATA + "explore-cgi/data/interim/cgi_ind_exp/pysster_fa/both.indel.sample.fa.gz"]
#chroms = []
#starts = []
#ends = []
#chroms, starts, ends = get_fa_info(fa_files)
#pysster_classifications = pd.DataFrame({'chrom' : list(chroms), 'pos_start' : list(starts), 'pos_end' : list(ends), 'labels' : list(labels), 'predictions' : list(predictions)}, columns = ['chrom','pos_start','pos_end','labels','predictions'])
#pysster_classifications.to_csv(DATA + 'explore-cgi/data/interim/pysster_classification_results/' + class_file_name + '.csv', sep = '\t', header = False, index = False)

utils.plot_roc(labels, predictions, output_folder+"roc.png")
utils.plot_prec_recall(labels, predictions, output_folder+"prec.png")
print(utils.get_performance_report(labels, predictions))

Image(output_folder+"roc.png")
Image(output_folder+"prec.png")

activations = model.get_max_activations(data, "all")
logos = model.visualize_all_kernels(activations, data, output_folder)
Image(output_folder+"motif_kernel_13.png")
Image(output_folder+"activations_kernel_13.png")
Image(output_folder+"position_kernel_13.png")
Image(output_folder+"data/alu.png")

utils.save_as_meme([logo[0] for logo in logos], output_folder+"motifs_seq.meme")
utils.save_as_meme([logo[1] for logo in logos], output_folder+"motifs_struct.meme")
model.plot_clustering(activations, output_folder+"clustering.png")
Image(output_folder+"clustering.png")

