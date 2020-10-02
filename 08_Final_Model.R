#' Title: Build the final version of PCAmodel from pre-processed data
#' Input: 
#' 1. `_PCs_rowNorm.rds`: An output from `05_PCA.R`, which is a list with the length 
#' of training data. Each element contains `rotation` and `variance`.
#' 2. `_PCclusters_hclust.rds`: An output from `06_Clustering.R`, which is a list 
#' of 7 elements.
#' 3. Other input parameters: `traningDatasets` and `note` about the model
#' 
#' Output: The final model named with suffix, `_PCAmodel.rds`
#' 
#' Process:
#' - 





library(PCAGenomicSignatures)
library(dplyr)

## Input parameters for PCAmodel_536
trainingDatasets <- "refinebioRseq"
note <- "536 refine.bio studies/ use top 20 PCs/ top 90% varying genes"

## Working directory
wd <- file.path("~/data2/PCAGenomicSignatureLibrary", 
                trainingDatasets, "PCAmodel_536")

## Load the training dataset information
dir <- system.file("extdata", package = "PCAGenomicSignatures")
studyMeta <- read.table(file.path(dir, "studyMeta.tsv"))
ind <- which(studyMeta$PCAmodel_536 == TRUE)
allStudies <- studyMeta$studyName[ind]

## Load PCA and Clustering results
trainingData_PCA <- readRDS(file.path(wd, paste0(trainingDatasets, "_PCs_rowNorm.rds")))  # PCA result
fname <- paste0(trainingDatasets, "_PCclusters.rds")
trainingData_PCclusters <- readRDS(file.path(wd, fname))  # Clustering result

## Variance Explained from PCA result
pca_summary <- list()
for (i in seq_along(trainingData_PCA)) {
    pca_summary[[i]] <- trainingData_PCA[[i]]$variance
    names(pca_summary)[i] <- names(trainingData_PCA)[i]
}

##### MeSH terms ###############################################################
dir <- system.file("extdata", package = "PCAGenomicSignatures")
x <- readRDS(file.path(dir, "MeSH_terms_1399refinebio.rds"))

# Update bagOfWords and MeSH_freq
MeSH_freq <- table(x$name) %>% sort(., decreasing = TRUE)  # freq. of each term
for (i in 1:nrow(x)) {x$bagOfWords[i] <- MeSH_freq[x$name[i]]}  # add freq. to the main table

# Order based on the score
x$score <- as.numeric(x$score)
x <- x[order(x$score, decreasing = TRUE),]

# 1398 studies in the given table
unique_id <- unique(x$identifier)

# Split MeSH term table to a list of each study, `all_MeSH`
all_MeSH <- vector("list", length(unique_id))
names(all_MeSH) <- unique_id
for (study in unique_id) {
    ind <- grepl(study, x$identifier)
    all_MeSH[[study]] <- x[ind, c("score", "identifier", "name", "explanation", "bagOfWords")]
}

# Subset only to the studies used in the model 
trainingData_MeSH <- all_MeSH[allStudies]


##### Build PCAGenomicSignatures object ########################################
PCAmodel <- PCAGenomicSignatures(assays = list(model = as.matrix(trainingData_PCclusters$avgLoading)))
metadata(PCAmodel) <- trainingData_PCclusters[c("cluster", "size", "k", "n")]
studies(PCAmodel) <- trainingData_PCclusters$studies
silhouetteWidth(PCAmodel) <- trainingData_PCclusters$sw
metadata(PCAmodel)$MeSH_freq <- MeSH_freq
trainingData(PCAmodel)$PCAsummary <- pca_summary
mesh(PCAmodel) <- trainingData_MeSH
updateNote(PCAmodel) <- note

# # GSEA
# source('~/data2/PCAGenomicSignatures/inst/script/build_gsea_DB.R')
# source('~/data2/PCAGenomicSignatures/inst/script/searchPathways.R')
# gsea_all <- searchPathways(PCAmodel, "~/data2/PCAGenomicSignatureLibrary/refinebioRseq/PCAmodel_536/gsea")
# gsea(PCAmodel) <- gsea_all


fname <- paste0(trainingDatasets, "_PCAmodel.rds")
saveRDS(PCAmodel, file.path(wd, fname))
