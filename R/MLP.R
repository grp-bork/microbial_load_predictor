#' Function to predict fecal microbial load based on species-level taxonomic profile of the human gut microbiome.
#' This function accepts the defaults outputs by mOTUs (v2.5 and v3), metaphlan (v3 and v4) profilers for shotgun metagenomes, or DADA2 with RDP classifier for 16S rRNA gene data.
#' @title MLP
#' @importFrom stats predict
#' @importFrom dplyr %>%
#' @importFrom here here
#' @importFrom readr read_rds
#' @importFrom vegan diversity
#' @param input Species-level (metagenome) or genus-level (16S rRNA) taxonomic profile of the gut microbiome.
#' @param profiler A profiler and its version used to generate the input specie-level taxonomic profile. "motus25", "motus3", "metaphlan3", "metaphlan4_mpa_vJan21_CHOCOPhlAnSGB_202103", "metaphlan4_mpa_vJun23_CHOCOPhlAnSGB_202307", or "rdp_train_set_16 (16S rRNA)". 
#' @param output Output format. "load": predicted microbial loads, or "qmp": taxonomic profile that takes into account the predicted loads. 

MLP <- function(input, profiler = "motus25", training_data = "metacardis", output = "load"){
  suppress_all <- function(expr) {
    suppressWarnings(suppressMessages(capture.output(expr)))
  }
  
  ## read prediction model
  ## MetaCardis model
  if(grepl("metacardis", training_data)){
    if(grepl("motus25", profiler)){
      model.path <- "data/metacardis/model.motus25.rds"
    }
    if(grepl("motus3", profiler)){
      model.path <- "data/metacardis/model.motus3.rds"
    }
    if(grepl("metaphlan3", profiler)){
      model.path <- "data/metacardis/model.metaphlan3.rds"
    }
    if(grepl("metaphlan4_mpa_vJan21_CHOCOPhlAnSGB_202103", profiler)){
      model.path <- "data/metacardis/model.metaphlan4.mpa_vJan21_CHOCOPhlAnSGB_202103.rds"
    }
    if(grepl("metaphlan4_mpa_vJun23_CHOCOPhlAnSGB_202307", profiler)){
      model.path <- "data/metacardis/model.metaphlan4.mpa_vJun23_CHOCOPhlAnSGB_202307.rds"
    }
  }
  
  ## GALAXY/MicrobLiver model  
  if(grepl("galaxy", training_data)){
    if(grepl("motus25", profiler)){
      model.path <- "data/galaxy/model.motus25.rds"
    }
    if(grepl("motus3", profiler)){
      model.path <- "data/galaxy/model.motus3.rds"
    }
    if(grepl("metaphlan3", profiler)){
      model.path <- "data/galaxy/model.metaphlan3.rds"
    }
    if(grepl("metaphlan4_mpa_vJan21_CHOCOPhlAnSGB_202103", profiler)){
      model.path <- "data/galaxy/model.metaphlan4.mpa_vJan21_CHOCOPhlAnSGB_202103.rds"
    }
    if(grepl("metaphlan4_mpa_vJun23_CHOCOPhlAnSGB_202307", profiler)){
      model.path <- "data/galaxy/model.metaphlan4.mpa_vJun23_CHOCOPhlAnSGB_202307.rds"
    }
  }
  
  ## 16S rRNA gene model
  if(grepl("rdp_train_set_16", profiler)){
    model.path <- "data/16S_rRNA/model.16S_rRNA.rds"
  }

  cat("Model Path:", model.path, "\n")
  model <- read_rds(model.path)

  ## change "unassigned" to "-1" (mOTUs v3.0)
  colnames(input) <- colnames(input) %>% str_replace("unassigned", "-1")

  ## add Shannon diversity
  input$`Shannon diversity` <- diversity(input)

  ## check colnames and add missing species if necessary
  d.tr <- model$trainingData
  d.tr <- d.tr[, -ncol(d.tr)]
  
  keep <- colnames(d.tr) %in% colnames(input)
  missing_sp <- colnames(d.tr)[!keep]
  for(i in missing_sp){
    input[[i]] <- 0
  }
  input <- input[, match(colnames(d.tr), colnames(input))]
  
  ## show the number of species detected in the input
  species_pro <- round(100 * sum(keep) / ncol(d.tr))
  sprintf(
    "\nℹ️ INFO: %d species were used in the selected model, and %d (%d%%) were found in the input file. Missing species have been supplemented.\n",
    ncol(d.tr), sum(keep), species_pro
  ) %>% message()
  
  ## predict microbial loads  
  input.min <- min(input[input != 0])/2
  temp <- suppress_all(
    pred <- predict(model, scale(log10(input + input.min)))
  )
  
  ## output is microbial load
  if(grepl("load", output)){
    res <- data.frame(load = 10 ^ pred)
    rownames(res) <- rownames(input)
    res <- res %>% rownames_to_column(var = "sample ID")
    return(res)
  }
  
  ## output is quantitative microbiome profile
  if(grepl("qmp", output)){
    keep <- !grepl("Shannon diversity", colnames(input))
    input <- input[, keep]
    q <- input * 10 ^ pred
    q <- q %>% rownames_to_column(var = "sample ID")
    return(q)
  }
}
