#' Title: Perform PCA and collect top 20 PCs
#' 
#' Input: 
#' 1. `cg` : An output from `04_Subset_Genes.R`, which is a list of common, top varying
#' genes from all the training datasets. File named with prefix `topGenes_`
#' 2. `in.dir`: Outputs from `02_Filtering.R`, which are the log-transformed count matrix of
#' samples in each study. Files named with suffix `_count.csv`.
#' 3. Others: `n` (the number of top PCs to collect) and `trainingDataSets` (the 
#' name of the training datasets)
#' 
#' Output: 
#' 1. `.rds` : a matrix containing all samples (with a subset of genes)
#' 2. `_SdMean.rds` : standard deviation and mean (with a subset of genes)
#' 3. A list with the length of studies in training data, names with suffix
#' `_PCs_rowNorm.rds`. Each element is a list of length 2 containing PCA results,
#' *rotation* and *variance*.
#' - rotation: a matrix with genes x top PCs 
#' - variance: a matrix with three rows (SD, Variance, Cumulative) x top PCs
#'  
#' Process:
#' 1. Subset each count matrix to the top varying, common genes
#' 2. Calculate mean and standard deviation from all samples and save in `_SdMean.rds`
#' 3. Row normalization
#' 4. Perform PCA: save top PCs and variance explained in `_PCs_rowNorm.rds`




## Input parameters for RAVmodel_536
n <- 20   # The number PCs to keep
cg <- readRDS("data/topGenes_13934.rds")
trainingDatasets <- "refinebioRseq" 

## Working directories
in.dir <- "~/data2/refinebio_processed/rna_seq_v2"
# for output
wd_data <- file.path("~/data2/GenomicSuperSignatureLibrary", trainingDatasets)  
if (!dir.exists(wd_data)) {dir.create(wd_data)}
wd <- file.path(wd_data, "RAVmodel_536")
if (!dir.exists(wd)) {dir.create(wd)}

## Load the training dataset information
dir <- system.file("extdata", package = "GenomicSuperSignature")
studyMeta <- read.table(file.path(dir, "studyMeta.tsv.gz"))
ind <- which(studyMeta$RAVmodel_536 == TRUE)
allStudies <- studyMeta$studyName[ind]

## An empty list for PCA results (rotation and variance)
trainingData_PCA <- vector("list", length(allStudies))
names(trainingData_PCA) <- allStudies

##### Calculate sd and mean across all the samples #############################
allSamples <- data.frame(matrix(NA, nrow = length(cg))) 
rownames(allSamples) <- cg
iter <- 0

for (study in allStudies) {
  dir.path <- file.path(in.dir, study)
  x <- data.table::fread(file.path(dir.path, paste0(study, "_count.csv")), 
                         showProgress = FALSE)
  x <- data.frame(x[,-1], row.names = x$V1) %>% as.matrix
  iter <- iter + 1
  allSamples <- cbind(allSamples, x[cg,])
}

allSamples <- allSamples[,-1]
s <- apply(allSamples, 1, sd)
m <- apply(allSamples, 1, mean)
SdMean <- list(sd = s, mean = m)
fname <- paste0(trainingDatasets, "_", iter, "study")

saveRDS(SdMean, file.path(wd, paste0(fname, "_SdMean.rds")))   # sd and mean
saveRDS(allSamples, file.path(wd, paste0(fname, ".rds")))   # matrix containing all samples
        
## Remove non-expressing genes in all samples (m == 0)
non_exp <- which(s == 0) %>% names
s <- s[!names(s) %in% non_exp]
m <- m[!names(m) %in% non_exp]

##### PCA ######################################################################
for (study in allStudies) {
  dir.path <- file.path(in.dir, study)
  x <- data.table::fread(file.path(dir.path, paste0(study, "_count.csv")), 
                         showProgress = FALSE)
  x <- data.frame(x[,-1], row.names = x$V1) %>% as.matrix
  x <- x[cg,,drop=FALSE]
  
  # Remove non-expressing genes in all samples (m == 0)
  x <- x[!rownames(x) %in% non_exp,]

  # This part will be used if any sample is removed upon filtering
  if (ncol(x) <= 20) {
    print(paste(study, "has only", ncol(x), "samples after filtering."))
    next
  }

  # Normalization
  x <- sweep(x, 1, m)
  x <- sweep(x, 1, s, "/")
  # x <- preprocessCore::normalize.quantiles(x)   # quantile normalization
  
  # PCA
  pca_res <- prcomp(t(x))   # x is a matrix with genes(row) x samples(column)
  trainingData_PCA[[study]]$rotation <- pca_res$rotation[,1:n]
  colnames(trainingData_PCA[[study]]$rotation) <- paste0(study, ".PC", c(1:n))
  eigs <- pca_res$sdev^2
  pca_summary <- rbind(SD = sqrt(eigs),
                       Variance = eigs/sum(eigs),
                       Cumulative = cumsum(eigs)/sum(eigs))
  trainingData_PCA[[study]]$variance <- pca_summary[,1:n]
  colnames(trainingData_PCA[[study]]$variance) <- paste0(study, ".PC", c(1:n))
  
  rm(x)
}

## Save
fname <- paste0(trainingDatasets, "_PCs_rowNorm.rds")
saveRDS(trainingData_PCA, file.path(wd, fname))
