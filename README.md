# Microbial load predictor (MLP)

R function to predict the fecal microbial load (total microbial cell count per gram or cell density) based on the taxonomic profile of the human gut microbiome. The prediction model was trained and constructed based on paired data of fecal metagenomes and fecal microbial load in the GALAXY and MetaCardis projects. Please also see the website for more information on this tool https://microbiome-tools.embl.de/mlp/.

- GALAXY dataset (Nishijima S et al, 2024, bioRxiv)  
[Fecal microbial load is a major determinant of gut microbiome variation and a confounder for disease associations](https://www.biorxiv.org/content/10.1101/2024.03.18.584290v1)

- MetaCardis dataset (Forslund, SK et al, 2021)  
[Combinatorial, additive and dose-dependent drug–microbiome associations](https://www.nature.com/articles/s41586-021-04177-9)

## Requirements
- R 4.3.1+  
- vegan
- tidyverse
- here

## Input file
Species-level taxonomic profiles prepared by the following taxonomic profilers (the default output) are supported. Example files from these profilers are available in the `test_data` folder.  

- mOTUs v2.5 (Milanese A et al., 2019)  
[Microbial abundance, activity and population genomic profiling with mOTUs2](https://www.nature.com/articles/s41467-019-08844-4)

- mOTUs v3.0 (Ruscheweyh HJ et al., 2022)  
[Cultivation-independent genomes greatly expand taxonomic-profiling capabilities of mOTUs across various environments](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-022-01410-z)

- MetaPhlAn3 (Beghini F et al., 2021) (mpa_v30_CHOCOPhlAn_201901)  
[Integrating taxonomic, functional, and strain-level profiling of diverse microbial communities with bioBakery 3](https://elifesciences.org/articles/65088)

- MetaPhlAn4 (Blanco-Míguez A et al., 2023) (mpa_vJan21_CHOCOPhlAnSGB_202103)  
[Extending and improving metagenomic taxonomic profiling with uncharacterized species using MetaPhlAn 4](https://www.nature.com/articles/s41587-023-01688-w)


## How to start
- On terminal
```
git clone https://git.embl.de/grp-bork/microbial-load-predictor.git
```

- On R (inside the downloaded folder)
```
devtools::install()
library("MLP")
```

## Predicting microbial load
```
load <- MLP(input, "motus2", "load")
```

## Transforming relative microbiome profile (RMP) to quantitative microbiome profile (QMP)
```
qmp <- MLP(input, "motus2", "qmp")
```
Quantitative (absolute) abundance = relative abundance * predicted microbial load

## Example code using test data
The test data comes from `Franzosa EA et al., 2018` including Crohn's disease and ulcerative colitis patients as well as control individuals.  
[Gut microbiome structure and metabolic activity in inflammatory bowel disease](https://www.nature.com/articles/s41564-018-0306-4)

```
library(tidyverse)

# read input file (mOTUs v2.5)
input <- read.delim("test_data/Franzosa_2018_IBD.motus25.tsv", header = T, row.names = 1, check.names = F) 

# transpose the data
input <- data.frame(t(input), check.names = F)

# predict microbial loads
load <- MLP(input, "motus2", "load")

# transform relative microbiome profile (RMP) to quantitative microbiome profile (QMP)
qmp <- MLP(input, "motus2", "qmp")

# plot predicted microbial loads (ggplot2 is required)
md <- read.delim("test_data/Franzosa_2018_IBD.metadata.tsv", header = T, row.names = 1, check.names = F)
df <- data.frame(md, load = load$load)
ggplot(df, aes(x = Disease, y = log10(load), fill = Disease)) +
  theme_bw() +
  geom_boxplot()
```
