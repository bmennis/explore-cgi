input_file <- "data/raw/test_set.vcf"
output_file <- "data/interim/test_set_output"

library(ggplot2)

d = read.delim(input_file, header=TRUE, sep='\t')
p = ggplot(data=d, aes(x=2)) + geom_histogram(binwidth=5)
ggsave(output_file, p)



