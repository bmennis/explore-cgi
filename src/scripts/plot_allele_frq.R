library(ggplot2)
          
d <- read.delim("data/raw/Kaviar_subset_allele_frq.csv", header=TRUE, sep='\t')
p <- ggplot(data=d) + geom_violin(aes(x=Sources, y=log(Allele.Frequency,10))) + theme_bw(base_size=18) + xlab("Sources") + ylab("Allele Frequencies")
ggsave('allele_plot.png', p)

