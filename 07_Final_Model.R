#' Title: Build the final version of RAVmodel from pre-processed data
#' Input: 
#' 1. `_PCs_rowNorm.rds`: An output from `05_PCA.R`, which is a list with the length 
#' of training data. Each element contains `rotation` and `variance`.
#' 2. `_PCclusters_hclust.rds` (or `_PCclusters.rds`): An output from 
#' `06_Clustering.R`, which is a list of 7 elements.
#' 3. Other input parameters: `trainingDatasets` and `note` about the model
#' 
#' Output: The final model named with suffix, `_RAVmodel_{annotGeneSets}.rds`
#' 
#' Process:
#' 1. Collect the following information from pre-processing:
#'    - Variance explained by PCs
#'    - MeSH terms for each study
#'    - GSEA on each PCcluster
#' 2. Combine all the information and built PCAGenomicSignatures object




library(GenomicSuperSignature)
library(dplyr)
d <- 2.25

## Input parameters for RAVmodel_536
trainingDatasets <- "refinebioRseq"
note <- "536 refine.bio studies/ use top 20 PCs/ top 90% varying genes/ GSEA with MSigDB C2.all"
annotGeneSets <- "C2"   # will be used for the output name
annot_database <- "MSigDB C2"   # will be added into the RAVmodel metadata

## Working directory
wd <- file.path("~/data2/GenomicSuperSignatureLibrary", 
                trainingDatasets, "RAVmodel_536")
## Output directory for different cluster number
if (d != 2.25) {
    out_dir <- paste0(wd, "_clNum", d)
    if (!dir.exists(out_dir)) {dir.create(out_dir)}
} else {out_dir <- wd}

## Load the training dataset information
dir <- system.file("extdata", package = "GenomicSuperSignature")
studyMeta <- read.table(file.path(dir, "studyMeta.tsv.gz"))
ind <- which(studyMeta$RAVmodel_536 == TRUE)
allStudies <- studyMeta$studyName[ind] # character vector with NCBI study accession

## Load PCA and Clustering results
trainingData_PCA <- readRDS(file.path(wd, paste0(trainingDatasets, "_PCs_rowNorm.rds")))  # PCA result
fname <- paste0(trainingDatasets, "_PCclusters.rds")
trainingData_PCclusters <- readRDS(file.path(out_dir, fname))  # Clustering result

##### Variance Explained from PCA result #######################################
pca_summary <- list()
for (i in seq_along(trainingData_PCA)) {
    pca_summary[[i]] <- trainingData_PCA[[i]]$variance
    names(pca_summary)[i] <- names(trainingData_PCA)[i]
}

##### MeSH terms ###############################################################
dir <- system.file("extdata", package = "GenomicSuperSignaturePaper")
load(file.path(dir, "MeSH_terms_1399refinebio.rda"))
x <- mesh_table

# Update bagOfWords and MeSH_freq
MeSH_freq <- table(x$name) %>% sort(., decreasing = TRUE)  # freq. of each term
for (i in 1:nrow(x)) {x$bagOfWords[i] <- MeSH_freq[x$name[i]]}  # add freq. to the main table

# 1398 studies in the given table
unique_id <- unique(x$identifier)

# Split MeSH term table to a list of each study, `all_MeSH`
# RAVmodel is currently doesn't have 'concept'.
all_MeSH <- vector("list", length(unique_id))
names(all_MeSH) <- unique_id
for (study in unique_id) {
    ind <- grepl(study, x$identifier)
    all_MeSH[[study]] <- x[ind, c("score", "identifier", "concept", "name", 
                                  "explanation", "bagOfWords")]
}

# Subset only to the studies used in the model 
trainingData_MeSH <- all_MeSH[allStudies]




##### Build PCAGenomicSignatures object ########################################
RAVmodel <- PCAGenomicSignatures(assays = list(RAVindex = as.matrix(trainingData_PCclusters$avgLoading)))
metadata(RAVmodel) <- trainingData_PCclusters[c("cluster", "size", "k", "n")]
names(metadata(RAVmodel)$size) <- paste0("RAV", seq_len(ncol(RAVmodel)))
geneSets(RAVmodel) <- annot_database
studies(RAVmodel) <- trainingData_PCclusters$studies
silhouetteWidth(RAVmodel) <- trainingData_PCclusters$sw
metadata(RAVmodel)$MeSH_freq <- MeSH_freq
trainingData(RAVmodel)$PCAsummary <- pca_summary[rownames(trainingData(RAVmodel))]
mesh(RAVmodel) <- trainingData_MeSH[rownames(trainingData(RAVmodel))]
updateNote(RAVmodel) <- note
metadata(RAVmodel)$version <- "1.1.0"

## Save pre-GSEA version
fname <- paste0(trainingDatasets, "_RAVmodel_", 
                format(Sys.Date(), format="%Y%m%d"), ".rds")
saveRDS(RAVmodel, file.path(out_dir, fname))


##### GSEA #####################################################################
gsea_script <- "inst/scripts/build_gsea_DB.R"
source(gsea_script)  # This is the processing script. Doule-check the details.

dir2 <- system.file("scripts", package = "GenomicSuperSignature")
searchPathways_func <- file.path(dir2, "searchPathways.R")
source(searchPathways_func)  # load the function

gsea_dir <- file.path(out_dir, paste0("gsea_", annotGeneSets))  # GSEA C2 DB is saved here
gsea_all <- searchPathways(RAVmodel, gsea_dir)  
gsea(RAVmodel) <- gsea_all

## Save post-GSEA version
fname <- paste0(trainingDatasets, "_RAVmodel_", annotGeneSets, 
                "_", format(Sys.Date(), format="%Y%m%d"),".rds")
saveRDS(RAVmodel, file.path(out_dir, fname))

