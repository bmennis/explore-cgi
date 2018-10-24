import os
from pysster.Data import Data
from pysster import utils
from IPython.display import Image
DATA = "/mnt/isilon/dbhi_bfx/perry/brian/"
#establish output directory
output_folder =  DATA + "explore-cgi/data/interim/cgi_ind_exp/pysster_output/model_test_run_10_17_18_2_feats/"
if not os.path.isdir(output_folder):
    os.makedirs(output_folder)

#load the pysster prediction model
model = utils.load_model("/mnt/isilon/dbhi_bfx/perry/brian/explore_cgi/data/interim/cgi_ind_exp/pysster_output/run_10_17_18_2_feats/model.pkl")

add_cgi_features = [DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/cgi.indel.unsample__microsat.out"]

add_both_features = [x.replace('cgi.', 'both.') for x in add_cgi_features]

indel_len_feat = [DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/cgi.indel.unsample__indel_length.out",
                  DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.unsample__indel_length.out"]

#load the dataset as data
data = Data([DATA + "explore-cgi/data/interim/cgi_ind_exp/pysster_fa/cgi.indel.unsample.fa.gz", DATA + "explore-cgi/data/interim/cgi_ind_exp/pysster_fa/both.indel.unsample.fa.gz"], ("ACGT", "XDI"))

for x, y in zip(add_cgi_features, add_both_features):
        features = [x,y]
        data.load_additional_data(features, is_categorical=True)

data.load_additional_data(indel_len_feat, is_categorical=False)

#run the model of pysster on all of the data set
predictions = model.predict(data, "all")
predictions

labels = data.get_labels("all")
labels

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

