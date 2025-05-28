<img src="assets/MLP_logo2.png" alt="logo" width="250">

## Microbial load predictor (MLP)

R function to predict the fecal microbial load (total microbial cell count per gram or cell density) based on the taxonomic profile of the human gut microbiome. The prediction model was trained and constructed based on paired data of fecal metagenomes and fecal microbial load in the GALAXY/MicrobLiver and MetaCardis projects. Please also see the website for more information on this tool https://microbiome-tools.embl.de/mlp/.

- GALAXY dataset (Nishijima S et al, 2024)  
[Fecal microbial load is a major determinant of gut microbiome variation and a confounder for disease associations](https://www.sciencedirect.com/science/article/pii/S0092867424012042)

- MetaCardis dataset (Forslund, SK et al, 2021)  
[Combinatorial, additive and dose-dependent drug–microbiome associations](https://www.nature.com/articles/s41586-021-04177-9)

## Requirements
- R 4.3.1+  
- vegan
- tidyverse

## Input taxonomic profile
Input files are species-level taxonomic profiles for shotgun metagenomes and genus-level profiles for 16S rRNA gene data prepared by the following taxonomic profilers (the default output). Example files from these profilers are available in the `test_data` folder.  

- **mOTUs v2.5** (Milanese A et al., 2019)  
[Microbial abundance, activity and population genomic profiling with mOTUs2](https://www.nature.com/articles/s41467-019-08844-4)

- **mOTUs v3.0** (Ruscheweyh HJ et al., 2022)  
[Cultivation-independent genomes greatly expand taxonomic-profiling capabilities of mOTUs across various environments](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-022-01410-z)

- **MetaPhlAn3** (Beghini F et al., 2021) (mpa_v30_CHOCOPhlAn_201901)  
[Integrating taxonomic, functional, and strain-level profiling of diverse microbial communities with bioBakery 3](https://elifesciences.org/articles/65088)

- **MetaPhlAn4** (Blanco-Míguez A et al., 2023) (mpa_vJan21_CHOCOPhlAnSGB_20210, mpa_vJun23_CHOCOPhlAnSGB_202307, mpa_vJun23_CHOCOPhlAnSGB_202403, or mpa_vJan25_CHOCOPhlAnSGB_202503)  
[Extending and improving metagenomic taxonomic profiling with uncharacterized species using MetaPhlAn 4](https://www.nature.com/articles/s41587-023-01688-w)

- **DADA2 with RDP** (Blanco-Míguez A et al., 2023) (rdp_train_set_16)  
[High-resolution sample inference from Illumina amplicon data](https://www.nature.com/articles/nmeth.3869)


## How to start
- On terminal
```shell
git clone https://git.embl.de/grp-bork/microbial_load_predictor.git
```

- On R (inside the downloaded folder)
```R
devtools::install()
library("MLP")
```

## Main Function
The `MLP` function predicts microbial load from taxonomic profiles. Please also see the [documentation page](https://microbiome-tools.embl.de/mlp/documentation) for details on the input file format and model difference. 

### Usage:
```r
out <- MLP(input, profiler, training_data, output)
```

### Parameters:
| Parameter      | Description |
|--------------|-------------|
| `input`      | The input taxonomic profiles|
| `profiler`   | The taxonomic profiler used |
| `training_data` | The training dataset used for the model|
| `output`     | The type of output|

#### Options for taxonomic profilers
For shotgun metagenomes (species-level taxonomic profile)
- "motus25"
- "motus3"
- "metaphlan3"
- "metaphlan4_mpa_vJan21_CHOCOPhlAnSGB_202103"
- "metaphlan4_mpa_vJun23_CHOCOPhlAnSGB_202307"
- "metaphlan4.mpa_vJun23_CHOCOPhlAnSGB_202403"
- "metaphlan4.mpa_vJan25_CHOCOPhlAnSGB_202503"

For 16S rRNA gene data (genus-level taxonomic profile)
- "rdp_train_set_16" (DADA2 with RDP)

#### Options for training data used in the model
- "galaxy"
- "metacardis"

#### Options for output
- "load"
- "qmp"

### Predicting microbial load
```R
load <- MLP(input_profile, "motus25", "metacardis", "load")
```

### Transforming relative microbiome profile (RMP) to quantitative microbiome profile (QMP)
```R
qmp <- MLP(input_profile, "motus25", "metacardis", "qmp")
```
Quantitative (absolute) abundance = relative abundance * predicted microbial load

### Example code using test data
The test data comes from `Franzosa EA et al., 2018` including Crohn's disease and ulcerative colitis patients as well as control individuals.  
[Gut microbiome structure and metabolic activity in inflammatory bowel disease](https://www.nature.com/articles/s41564-018-0306-4)

```R
library(tidyverse)

# read input file (mOTUs v2.5)
input <- read.delim("test_data/Franzosa_2018_IBD.motus25.tsv", header = T, row.names = 1, check.names = F) 

# transpose the data
input <- data.frame(t(input), check.names = F)

# predict microbial loads with the MetaCardis model
load <- MLP(input, "motus25", "metacardis", "load")

# transform relative microbiome profile (RMP) to quantitative microbiome profile (QMP)
qmp <- MLP(input, "motus25", "metacardis", "qmp")

# plot predicted microbial loads (ggplot2 is required)
md <- read.delim("test_data/Franzosa_2018_IBD.metadata.tsv", header = T, row.names = 1, check.names = F)
df <- data.frame(md, load = load$load)
ggplot(df, aes(x = Disease, y = log10(load), fill = Disease)) +
  theme_bw() +
  geom_boxplot()
```
