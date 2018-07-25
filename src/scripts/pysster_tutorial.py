import os
from time import time
from IPython.display import Image
from pysster.Data import Data
from pysster.Grid_Search import Grid_Search
from pysster import utils

###Establish output directory
output_folder = "explore_cgi/data/interim/kaviar_pysster_additional_features/"
if not os.path.isdir(output_folder):
    os.makedirs(output_folder)


data = Data(["/home/ennisb/me/explore-cgi/data/interim/kaviar_fa_gz/short.cgi.fa.gz", "/home/ennisb/me/explore-cgi/data/interim/kaviar_fa_gz/short.both.fa.gz"], ("ACGT"))

add_cgi_features = ["/home/ennisb/me/explore-cgi/data/interim/additional_features/1kg.cgi.out","/home/ennisb/me/explore-cgi/data/interim/additional_features/20120824_combined_mask.cgi.out","/home/ennisb/me/explore-cgi/data/interim/additional_features/blackTerry.cgi.out","/home/ennisb/me/explore-cgi/data/interim/additional_features/dgv.cgi.out","/home/ennisb/me/explore-cgi/data/interim/additional_features/dgv.short.cgi.out","/home/ennisb/me/explore-cgi/data/interim/additional_features/GRCh37GenomicSuperDup.sorted.cgi.out","/home/ennisb/me/explore-cgi/data/interim/additional_features/hg19.blacklist.cgi.out","/home/ennisb/me/explore-cgi/data/interim/additional_features/rmsk.cgi.out","/home/ennisb/me/explore-cgi/data/interim/additional_features/simpleRepeat.cgi.out"]

add_both_features = ["/home/ennisb/me/explore-cgi/data/interim/additional_features/1kg.both.out","/home/ennisb/me/explore-cgi/data/interim/additional_features/20120824_combined_mask.both.out","/home/ennisb/me/explore-cgi/data/interim/additional_features/blackTerry.both.out","/home/ennisb/me/explore-cgi/data/interim/additional_features/dgv.both.out","/home/ennisb/me/explore-cgi/data/interim/additional_features/dgv.short.both.out","/home/ennisb/me/explore-cgi/data/interim/additional_features/GRCh37GenomicSuperDup.sorted.both.out","/home/ennisb/me/explore-cgi/data/interim/additional_features/hg19.blacklist.both.out","/home/ennisb/me/explore-cgi/data/interim/additional_features/rmsk.both.out","/home/ennisb/me/explore-cgi/data/interim/additional_features/simpleRepeat.both.out"]

for x, y in zip(add_cgi_features, add_both_features):
	features = [x,y]
	print(features)
	data.load_additional_data(features, is_categorical=True) 
	

print(data.get_summary())

data.train_val_test_split(portion_train=0.6, portion_val=0.2, seed=3)
print(data.get_summary())

###Model Training
params = {"conv_num": [2, 3], "kernel_num": [100], "kernel_len": [20], "dropout_input": [0.1, 0.4]}
searcher = Grid_Search(params)
start = time()
model, summary = searcher.train(data, pr_auc = True,  verbose=False)
stop = time()
print("time in minutes: {}".format((stop-start)/60))

print(summary)

###Perfomance evaluation
predictions = model.predict(data, "test")
predictions

labels = data.get_labels("test")
labels

utils.plot_roc(labels, predictions, output_folder+"roc.png")
utils.plot_prec_recall(labels, predictions, output_folder+"prec.png")
print(utils.get_performance_report(labels, predictions))

Image(output_folder+"roc.png")
Image(output_folder+"prec.png")

###Motif Visualization
activations = model.get_max_activations(data, "test")
logos = model.visualize_all_kernels(activations, data, output_folder)
Image(output_folder+"motif_kernel_13.png")
Image(output_folder+"activations_kernel_13.png")
Image(output_folder+"position_kernel_13.png")
Image("data/alu.png")

utils.save_as_meme([logo[0] for logo in logos], output_folder+"motifs_seq.meme")
utils.save_as_meme([logo[1] for logo in logos], output_folder+"motifs_struct.meme")
model.plot_clustering(activations, output_folder+"clustering.png")
Image(output_folder+"clustering.png")

utils.save_data(data, output_folder+"data.pkl")
utils.save_model(model, output_folder+"model.pkl")
