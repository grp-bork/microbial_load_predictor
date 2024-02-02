#' Function to predict fecal microbial load based on species-level taxonomic profile of the human gut microbiome.
#' This function accepts the defaults outputs by mOTUs (v2.5 and v3) and metaphlan (v3 and v4) profilers.
#' @title MLpredictor
#' @importFrom stats predict
#' @importFrom dplyr %>%
#' @importFrom here here
#' @importFrom readr read_rds
#' @importFrom vegan diversity
#' @param input Species-level taxonomic profile of the gut microbiome.
#' @param profiler A profiler and its version used to generate the specie-level taxonomic profile. "motus2", "motus3", "metaphlan3", or "metaphlan4". 
#' @param output Output format. "load": predicted microbial loads, or "qmp": taxonomic profile that takes into account the predicted loads. 

MLpredictor <- function(input, profiler = "motus2", output = "load"){
  
  ## read prediction model
  if(grepl("motus2", profiler)){
    model.path <- "data/model.motus25.rds" %>% here()
  }
  if(grepl("motus3", profiler)){
    model.path <- "data/model.motus3.rds" %>% here()
  }
  if(grepl("metaphlan3", profiler)){
    model.path <- "data/model.metaphlan3.rds" %>% here()
  }
  if(grepl("metaphlan4", profiler)){
    model.path <- "data/model.metaphlan4.rds" %>% here()
  }
  model <- read_rds(model.path)
  
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
  
  ## predict microbial loads  
  input.min <- min(input[input != 0])/2
  pred <- predict(model, scale(log10(input + input.min)))
  
  ## output is microbial load
  if(grepl("load", output)){
    res <- data.frame(load = 10 ^ pred)
    rownames(res) <- rownames(input)
    return(res)
  }
  
  ## output is quantitative microbiome profile
  if(grepl("qmp", output)){
    keep <- !grepl("Shannon diversity", colnames(input))
    input <- input[, keep]
    q <- input * 10 ^ pred
    return(q)
  }
}
