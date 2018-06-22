"""Scripts for plots."""

# plot af from kaivar on kaviar subsets
rule plot_kaviar_af_kaviar:
    input:  DATA + 'raw/Kaviar_subset_allele_frq.csv'
    output: FILES + 'plots/kavar_af.pdf'
    shell: 'Rscript {SCRITPS}plot_allele_frq.R {input} {output}'
