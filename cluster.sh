snakemake -c  "qsub -cwd -V -l h_vmem={cluster.h_vmem} -l mem_free={cluster.mem_free} -l m_mem_free={cluster.m_mem_free} -pe smp {threads}" \
         -j 500 --latency-wait 86400 --greediness 0.8 \
         --cluster-config configs/cluster.yml \
         -s src/rules/sf_ann_target.py mk_target_samp_add_feat mk_target_unsamp_add_feat

 
