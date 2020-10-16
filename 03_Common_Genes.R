#' Title: Collect all the genes from each study/experiment
#' Input: Count matrices from training datasets 
#' - input directory (`in.dir`): tximport-ed files are saved here. Same as the `out.dir`
#' from `01_Import.R` script.
#' 
#' Output: A list contining a list of genes from each study. The output file (named 
#' with the prefix `topGenesInTrainingData_`) is saved in the `data` folder within 
#' the directory this script is. Genes are listed in a decreasing order based on 
#' its standard deviation within each study. This level of variation is used later 
#' to subset genes for modeling buildilng. 
#' 
#' Process:
#' 1. Calculate standard deviation of each gene's expression at the study level
#' 2. Order genes from most varying to least varying
#' 3. Save ordered list of genes of each study in one object, `topGenesInTrainingData`




## Working directories
in.dir <- "/nobackup/16tb_b/GenomicSignatures/refinebio/rna_seq_v2"  
  
## Load the training dataset information
dir <- system.file("extdata", package = "PCAGenomicSignatures")
studyMeta <- read.table(file.path(dir, "studyMeta.tsv"))
ind <- which(studyMeta$PCAmodel_536 == TRUE)
studies <- studyMeta$studyName[ind]

## List containing the ordered list of genes from each study
topGenesInTrainingData <- list()  
for (study in studies) {
  dir.path <- file.path(in.dir, study)
  x <- data.table::fread(file.path(dir.path, paste0(study, "_count.csv")), 
                         showProgress = FALSE)
  
  # order by standard deviation
  y <- apply(x[,-1], 1, sd) 
  names(y) <- x$V1
  y <- y[order(y, decreasing = TRUE)]

  topGenesInTrainingData[[study]] <- y
  rm(x)
}

## Save
len <- length(topGenesInTrainingData)
fname <- paste0("topGenesInTrainingData_", len, "_", Sys.Date(), ".rds")
saveRDS(topGenesInTrainingData, file.path("data", fname))
