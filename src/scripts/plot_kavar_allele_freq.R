library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
dat_file = args[1]
plot_file = args[2

d <- read.delim(dat_file, header=TRUE, sep='\t')
p <- ggplot(data=d) + geom_violin(aes(x=Sources, y=log(Allele.Frequency,10))) + theme_bw(base_size=18) + xlab("Sources") + ylab("Allele Frequencies")
ggsave(plot_file, p)

