#' Title: Pre-process RNA-seq data from refine.bio 
#' Input: `tximpot` outputs (`.rds` format) of studies with > 20 downloaded samples  
#' - input directory (`in.dir`): tximport-ed files are saved here. Same as the `out.dir`
#' from `01_Import.R` script.
#' 
#' Output: One gene expression matrix for each study (= experiment) with the 
#' extension, `_count.csv`. Only log2(x+1) transformation is done. This file is
#' saved in the same directory as the input `.rds` file.
#' 
#' Process:
#' 1. Use 1,399 studies with > 20 samples succesfully downloaded and imported
#' 2. log2 transformation of the count matrix
#' 3. No filtering on any genes/samples
#' 
#' Note:
#' Some models has additional filtering criteria. For example, RAVmodel_536 further
#' excludes studies potentially from 'Single-Cell Analysis'. The final list of 
#' training datasets for each model is availble in the package through `extdata/studyMeta.tsv`.
#' Specifically, excluding 'Single-Cell Analysis' studies process is described 
#' here: https://shbrief.github.io/GenomicSuperSignaturePaper/Methods/training_datasets/available_samples.html#pcamodel_536.





library(dplyr)

## Working directories
in.dir <- "/nobackup/16tb_b/GenomicSignatures/refinebio/rna_seq_v2"  

## Load the training dataset information
dir <- system.file("extdata", package = "GenomicSuperSignature")
studyMeta <- read.table(file.path(dir, "studyMeta.tsv"))

## 1,399 studies that are successfully imported and have > 20 samples
ind <- which(studyMeta$imported == TRUE)
allStudies <- studyMeta$studyName[ind]

## Log2 transformation
for (study in allStudies) {
  dir.path <- file.path(in.dir, study)
  file.path <- file.path(dir.path, paste0(study, ".rds"))
  countMatrix <- readRDS(file.path)
  
  # log2(TPM+1) transformation
  count <- countMatrix$counts 
  count <- log2(count + 1)

  # Save the results for each study in `_count.csv`
  write.table(count, file.path(dir.path, paste0(study, "_count.csv")))

  # standard output to check the progress
  filteringDone <- paste("Successful Filtering done :", study)
  print(filteringDone)
}

Sys.time()
