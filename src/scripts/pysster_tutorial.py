import os
from time import time
from IPython.display import Image
from pysster.Data import Data
from pysster.Grid_Search import Grid_Search
from pysster import utils
DATA = '/mnt/isilon/cbmi/variome/perry/brian/'
###Establish output directory
output_folder = DATA + "explore_cgi/data/interim/kaviar_pysster_add_feat_indel/"
if not os.path.isdir(output_folder):
    os.makedirs(output_folder)


data = Data([DATA + "explore-cgi/data/interim/kaviar_fa_gz/full.cgi.indel.fa.gz", DATA + "explore-cgi/data/interim/kaviar_fa_gz/full.both.indel.fa.gz"], ("ACGT"))

add_cgi_features = [DATA + "explore-cgi/data/interim/additional_features/1kg.cgi.indel.out",
                    DATA + "explore-cgi/data/interim/additional_features/20120824_combined_mask.cgi.indel.out",
                    DATA + "explore-cgi/data/interim/additional_features/blackTerry.cgi.indel.out",
                    DATA + "explore-cgi/data/interim/additional_features/dgv.cgi.indel.out",
                    DATA + "explore-cgi/data/interim/additional_features/dgv.short.cgi.indel.out",
                    DATA + "explore-cgi/data/interim/additional_features/GRCh37GenomicSuperDup.sorted.cgi.indel.out",
                    DATA + "explore-cgi/data/interim/additional_features/hg19.blacklist.cgi.indel.out",
                    DATA + "explore-cgi/data/interim/additional_features/rmsk.cgi.indel.out",
                    DATA + "explore-cgi/data/interim/additional_features/simpleRepeat.cgi.indel.out",
                    DATA + "explore-cgi/data/interim/alignability/wgEncodeCrgMapabilityAlign100mer.cgi.indel.out",
                    DATA + "explore-cgi/data/interim/alignability/wgEncodeCrgMapabilityAlign24mer.cgi.indel.out",
                    DATA + "explore-cgi/data/interim/alignability/wgEncodeCrgMapabilityAlign36mer.cgi.indel.out",
                    DATA + "explore-cgi/data/interim/alignability/wgEncodeCrgMapabilityAlign40mer.cgi.indel.out",
                    DATA + "explore-cgi/data/interim/alignability/wgEncodeCrgMapabilityAlign50mer.cgi.indel.out",
                    DATA + "explore-cgi/data/interim/alignability/wgEncodeCrgMapabilityAlign75mer.cgi.indel.out",
                    DATA + "explore-cgi/data/interim/excludable/consensusBlacklist.cgi.indel.out",
                    DATA + "explore-cgi/data/interim/excludable/dukeExcludeRegions.cgi.indel.out",
                    DATA + "explore-cgi/data/interim/uniqueness/wgEncodeDukeMapabilityUniqueness20bp.cgi.indel.out",
                    DATA + "explore-cgi/data/interim/uniqueness/wgEncodeDukeMapabilityUniqueness35bp.cgi.indel.out"]

add_both_features = [x.replace('.cgi.', '.both.') for x in add_cgi_features]

indel_len_feat = [DATA + "explore-cgi/data/interim/additional_features/indel_length.cgi.indel.out",
                  DATA + "explore-cgi/data/interim/additional_features/indel_length.both.indel.out"]

for x, y in zip(add_cgi_features, add_both_features):
	features = [x,y]
	print(features)
	data.load_additional_data(features, is_categorical=True) 
	
data.load_additional_data(indel_len_feat, is_categorical=False)

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

##Perfomance evaluation
predictions = model.predict(data, "test")
predictions

labels = data.get_labels("test")
labels

utils.plot_roc(labels, predictions, output_folder+"roc.png")
utils.plot_prec_recall(labels, predictions, output_folder+"prec.png")
print(utils.get_performance_report(labels, predictions))

Image(output_folder+"roc.png")
Image(output_folder+"prec.png")


activations = model.get_max_activations(data, "test")
logos = model.visualize_all_kernels(activations, data, output_folder)
Image(output_folder+"motif_kernel_13.png")
Image(output_folder+"activations_kernel_13.png")
Image(output_folder+"position_kernel_13.png")
Image(output_folder+"data/alu.png")

#utils.save_as_meme([logo[0] for logo in logos], output_folder+"motifs_seq.meme")
#utils.save_as_meme([logo[1] for logo in logos], output_folder+"motifs_struct.meme")
model.plot_clustering(activations, output_folder+"clustering.png")
Image(output_folder+"clustering.png")

utils.save_data(data, output_folder+"data.pkl")
utils.save_model(model, output_folder+"model.pkl")
