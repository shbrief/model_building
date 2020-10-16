#' Title: Subset genes with different criteria
#' Input: List of ordered gene lists (based on standard deviation) from each study.
#' Files start with `topGenesInTrainingData_` prefix, followed by the number of 
#' studies and the date it's created.
#' - topGenesInTrainingData
#' - cutoff
#' 
#' Output: A list of the common genes met the filtering criteria. The output file 
#' (named with the prefix `topGenes_`) is saved in the `data` folder within the 
#' directory this script is.




## Choose the training dataset
topGenesInTrainingData <- readRDS("data/topGenesInTrainingData_536_2020-08-07.rds")

##### Top 90% varying genes ####################################################
cutoff <- 0.9  # Adjust this if you want different cutoff
topGenes <- c()

k <- length(topGenesInTrainingData[[1]])
k <- round(k*cutoff)
topGenes <- names(topGenesInTrainingData[[1]])[1:k]

for (i in 2:length(topGenesInTrainingData)) {
  l <- length(topGenesInTrainingData[[i]])
  l <- round(l*cutoff)
  ls <- list(topGenes, names(topGenesInTrainingData[[i]])[1:l])
  topGenes <- Reduce(intersect, ls)
}

saveRDS(topGenes, file.path("data", 
                            paste0("topGenes_", length(topGenes),".rds")))


##### Subset with priors #######################################################
# category <- "C2"
# m <- msigdbr::msigdbr(species = "Homo sapiens", category = category)
# cgPrior <- unique(m$gene_symbol)
# cgTrainingData <- Reduce(intersect, lapply(topGenesInTrainingData, names))
# cg <- intersect(cgPrior, cgTrainingData)
# saveRDS(cg, file.path("data", 
#                       paste0("topGenes_MSigDB", category, "_", length(cg),".rds")))

