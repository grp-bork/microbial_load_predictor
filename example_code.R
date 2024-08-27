## example script using test_data
library(tidyverse)
library(ggpubr)
source("R/MLP.R")

# read input file in the test_data (mOTUs v2.5)
input <- read.delim("test_data/Franzosa_2018_IBD.motus25.tsv", header = T, row.names = 1, check.names = F) 

# transpose the data
input <- data.frame(t(input), check.names = F)

# predict microbial loads
load <- MLP(input, "motus25", "metacardis", "load")
write_tsv(load, file = "~/Desktop/temp.tsv")

# transform relative microbiome profile (RMP) to quantitative microbiome profile (QMP)
qmp <- MLP(input, "motus25", "metacardis", "qmp")

# plot predicted microbial loads using ggplot2
md <- read.delim("test_data/Franzosa_2018_IBD.metadata.tsv", header = T, row.names = 1, check.names = F)
df <- data.frame(md, load = load$load)
comp <- combn(unique(df$Disease), 2, simplify = F)

ggplot(df, aes(x = Disease, y = log10(load), fill = Disease)) +
  theme_bw() +
  geom_boxplot() +
  ggpubr::stat_compare_means(comparisons = comp)
