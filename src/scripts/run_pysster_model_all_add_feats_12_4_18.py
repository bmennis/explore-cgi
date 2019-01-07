import os
from pysster.Data import Data
from pysster import utils
from IPython.display import Image
DATA = "/mnt/isilon/dbhi_bfx/perry/brian/"
#establish output directory
output_folder =  DATA + "explore-cgi/data/interim/cgi_ind_exp/pysster_output/model_run_12_5_18_tar_cgi_kav_both/"
if not os.path.isdir(output_folder):
    os.makedirs(output_folder)

#load the pysster prediction model
model = utils.load_model("/mnt/isilon/dbhi_bfx/perry/brian/explore_cgi/data/interim/cgi_ind_exp/pysster_output/train_run_10_18_18_all_add_feats_back/model.pkl")

add_cgi_features = [DATA + "explore-cgi/data/interim/target_exp/add_feats/cgi.indel.microsat.out",
                    DATA + "explore-cgi/data/interim/target_exp/add_feats/cgi.indel.lowmappabilityall.out",
                    DATA + "explore-cgi/data/interim/target_exp/add_feats/cgi.indel.notinlowmappabilityall.out",
                    DATA + "explore-cgi/data/interim/target_exp/add_feats/cgi.indel.siren_similarRegions_dist1.out",
                    DATA + "explore-cgi/data/interim/target_exp/add_feats/cgi.indel.segdupall.out",
                    DATA + "explore-cgi/data/interim/target_exp/add_feats/cgi.indel.notinsegdupall.out",
                    DATA + "explore-cgi/data/interim/target_exp/add_feats/cgi.indel.notinrefseq_union_cds.sort.out",
                    DATA + "explore-cgi/data/interim/target_exp/add_feats/cgi.indel.BadPromoters_gb-2013-14-5-r51-s1.out",
                    DATA + "explore-cgi/data/interim/target_exp/add_feats/cgi.indel.human_g1k_v37_l100_gclt30orgt55_slop50.out",
                    DATA + "explore-cgi/data/interim/target_exp/add_feats/cgi.indel.human_g1k_v37_l100_gc30to55_slop50.out",
                    DATA + "explore-cgi/data/interim/target_exp/add_feats/cgi.indel.hg19_self_chain_split.sort.out",
                    DATA + "explore-cgi/data/interim/target_exp/add_feats/cgi.indel.hg19_self_chain_split_both.out",
                    DATA + "explore-cgi/data/interim/target_exp/add_feats/cgi.indel.notinAllRepeats_gt95percidentity_slop5.out",
                    DATA + "explore-cgi/data/interim/target_exp/add_feats/cgi.indel.AllRepeats_gt95percidentity_slop5.out",
                    DATA + "explore-cgi/data/interim/target_exp/add_feats/cgi.indel.AllRepeats_lt51bp_gt95identity_merged.out",
                    DATA + "explore-cgi/data/interim/target_exp/add_feats/cgi.indel.1kg.out",
                    DATA + "explore-cgi/data/interim/target_exp/add_feats/cgi.indel.20120824_combined_mask.out",
                    DATA + "explore-cgi/data/interim/target_exp/add_feats/cgi.indel.blackTerry.out",
                    DATA + "explore-cgi/data/interim/target_exp/add_feats/cgi.indel.dgv.out",
                    DATA + "explore-cgi/data/interim/target_exp/add_feats/cgi.indel.dgv.short.out",
                    DATA + "explore-cgi/data/interim/target_exp/add_feats/cgi.indel.GRCh37GenomicSuperDup.sorted.out",
#                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/cgi.indel.sample__hg19.blacklist.out",
                    DATA + "explore-cgi/data/interim/target_exp/add_feats/cgi.indel.rmsk.out",
                    DATA + "explore-cgi/data/interim/target_exp/add_feats/cgi.indel.simpleRepeat.out",
#                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/cgi.indel.sample__wgEncodeCrgMapabilityAlign100mer.out",
#                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/cgi.indel.sample__wgEncodeCrgMapabilityAlign24mer.out",
#                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/cgi.indel.sample__wgEncodeCrgMapabilityAlign36mer.out",
#                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/cgi.indel.sample__wgEncodeCrgMapabilityAlign40mer.out",
#                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/cgi.indel.sample__wgEncodeCrgMapabilityAlign50mer.out",
#                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/cgi.indel.sample__wgEncodeCrgMapabilityAlign75mer.out",
#                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/cgi.indel.sample__consensusBlacklist.out",
#                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/cgi.indel.sample__dukeExcludeRegions.out",
#                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/cgi.indel.sample__wgEncodeDukeMapabilityUniqueness20bp.out",
#                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/cgi.indel.sample__wgEncodeDukeMapabilityUniqueness35bp.out",
                    DATA + "explore-cgi/data/interim/target_exp/add_feats/cgi.indel.SimpleRepeat_homopolymer_6to10.out",
                    DATA + "explore-cgi/data/interim/target_exp/add_feats/cgi.indel.SimpleRepeat_homopolymer_gt10.out",
                    DATA + "explore-cgi/data/interim/target_exp/add_feats/cgi.indel.SimpleRepeat_imperfecthomopolgt10_slop5.out",
                    DATA + "explore-cgi/data/interim/target_exp/add_feats/cgi.indel.GSE63874_Na_K_minus_hits_intersect.out",
                    DATA + "explore-cgi/data/interim/target_exp/add_feats/cgi.indel.GSE63874_Na_PDS_minus_hits_intersect.out",
                    DATA + "explore-cgi/data/interim/target_exp/add_feats/cgi.indel.Na_K_plus_hits_intersect.out",
                    DATA + "explore-cgi/data/interim/target_exp/add_feats/cgi.indel.Na_PDS_plus_hits_intersect.out",
                    DATA + "explore-cgi/data/interim/target_exp/add_feats/cgi.indel.hg19.2014.regions.out",
                    DATA + "explore-cgi/data/interim/target_exp/add_feats/cgi.indel.interrupted_repeats.out",
                    DATA + "explore-cgi/data/interim/target_exp/add_feats/cgi.indel.cpg_island_ext.out"]

add_both_features = [DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.sample__microsat.out",
                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.sample__lowmappabilityall.out",
                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.sample__notinlowmappabilityall.out",
                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.sample__siren_similarRegions_dist1.out",
                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.sample__segdupall.out",
                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.sample__notinsegdupall.out",
                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.sample__notinrefseq_union_cds.sort.out",
                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.sample__BadPromoters_gb-2013-14-5-r51-s1.out",
                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.sample__human_g1k_v37_l100_gclt30orgt55_slop50.out",
                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.sample__human_g1k_v37_l100_gc30to55_slop50.out",
                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.sample__hg19_self_chain_split.sort.out",
                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.sample__hg19_self_chain_split_both.out",
                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.sample__notinAllRepeats_gt95percidentity_slop5.out",
                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.sample__AllRepeats_gt95percidentity_slop5.out",
                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.sample__AllRepeats_lt51bp_gt95identity_merged.out",
                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.sample__1kg.out",
                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.sample__20120824_combined_mask.out",
                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.sample__blackTerry.out",
                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.sample__dgv.out",
                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.sample__dgv.short.out",
                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.sample__GRCh37GenomicSuperDup.sorted.out",
#                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/cgi.indel.unsample__hg19.blacklist.out",
                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.sample__rmsk.out",
                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.sample__simpleRepeat.out",
#                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/cgi.indel.unsample__wgEncodeCrgMapabilityAlign100mer.out",
#                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/cgi.indel.unsample__wgEncodeCrgMapabilityAlign24mer.out",
#                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/cgi.indel.unsample__wgEncodeCrgMapabilityAlign36mer.out",
#                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/cgi.indel.unsample__wgEncodeCrgMapabilityAlign40mer.out",
#                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/cgi.indel.unsample__wgEncodeCrgMapabilityAlign50mer.out",
#                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/cgi.indel.unsample__wgEncodeCrgMapabilityAlign75mer.out",
#                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/cgi.indel.unsample__consensusBlacklist.out",
#                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/cgi.indel.unsample__dukeExcludeRegions.out",
#                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/cgi.indel.unsample__wgEncodeDukeMapabilityUniqueness20bp.out",
#                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/cgi.indel.unsample__wgEncodeDukeMapabilityUniqueness35bp.out",
                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.sample__SimpleRepeat_homopolymer_6to10.out",
                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.sample__SimpleRepeat_homopolymer_gt10.out",
                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.sample__SimpleRepeat_imperfecthomopolgt10_slop5.out",
                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.sample__GSE63874_Na_K_minus_hits_intersect.out",
                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.sample__GSE63874_Na_PDS_minus_hits_intersect.out",
                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.sample__Na_K_plus_hits_intersect.out",
                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.sample__Na_PDS_plus_hits_intersect.out",
                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.sample__hg19.2014.regions.out",
                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.sample__interrupted_repeats.out",
                    DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.sample__cpg_island_ext.out"]

#add_both_features = [x.replace('cgi.indel.sample', 'both.indel.sample') for x in add_both_features]

indel_len_feat = [DATA + "explore-cgi/data/interim/target_exp/add_feats/cgi.indel.indel_length.out",
                  DATA + "explore-cgi/data/interim/cgi_ind_exp/add_feat/both.indel.sample__indel_length.out"]

#load the dataset as data
data = Data([DATA + "explore-cgi/data/interim/target_exp/target_pysster_fa/cgi.indel.fa.gz", DATA + "explore-cgi/data/interim/cgi_ind_exp/pysster_fa/both.indel.sample.fa.gz"], ("ACGT", "XDI"))

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

