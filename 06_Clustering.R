#' Title: Hierarchical Clustering of PCs from training data
#' Input: A list with the length of training datasets, named with the suffix 
#' `_PCs_rowNorm.rd `. Each element comtains `rotation` and `variance`.
#' Output: A list of 7 elements, including avgLoading. Output files named with 
#' the suffix `_PCclusters.rds`.
#' 
#' Process:
#' 1. Combine all the PCs
#' 2. Calculate distance matrix based on Spearman correlation
#' 3. Hierarchical clustering and cut the tree at `k = round(#ofPCs/d, 0)`
#' 4. Create avgLoading using `GenomicSuperSignature::buildAvgLoading`
#' 5. Calculate Silhouette width of each cluster





library(GenomicSuperSignature)
library(factoextra)
trainingDatasets <- "refinebioRseq"
d <- 4  # for cluster number (line 56)

## Working directory
wd <- file.path("~/data2/GenomicSuperSignatureLibrary", 
                trainingDatasets, "RAVmodel_536")

## Output directory for different cluster number
if (d != 2.25) {
    out_dir <- paste0(wd, "_clNum", d)
    if (!dir.exists(out_dir)) {dir.create(out_dir)}
} else {out_dir <- wd}

## Import input data
fname <- paste0(trainingDatasets, "_PCs_rowNorm.rds")
trainingData_PCA <- readRDS(file.path(wd, fname))
print("PCs.rds file is loaded")

## Combine all PCs
allZ_list <- lapply(trainingData_PCA, function(x) x$rotation)
allZ <- Reduce(cbind, allZ_list)
all <- t(allZ)   # a matrix of PCs (row) x genes (column)
print(paste("Dimension of allZ is", dim(allZ)))   # 13,934 genes x 10,720 PCs
saveRDS(allZ, file.path(out_dir, "allZ.rds"))


##### Hierarchical Clustering ##################################################
## Calculate distance
start <- Sys.time()
res.dist <- factoextra::get_dist(all, method = "spearman")
end <- Sys.time()
t <- end - start
print(t)  
saveRDS(res.dist, file.path(out_dir, "res_dist.rds"))

## Cut the tree
k <- round(nrow(all)/d, 0)
start <- Sys.time()
res.hcut <- factoextra::hcut(res.dist, k = k, hc_func = "hclust", 
                             hc_method = "ward.D", hc_metric = "spearman")
end <- Sys.time()
t2 <- end - start
print(t2)   
saveRDS(res.hcut, file.path(out_dir, "res_hcut.rds"))


##### Build avgLoading #########################################################
trainingData_PCclusters <- buildAvgLoading(allZ, k, cluster = res.hcut$cluster)

## Silhouette Width
cl <- trainingData_PCclusters$cluster
silh_res <- cluster::silhouette(cl, res.dist)
cl_silh_width <- summary(silh_res)$clus.avg.widths
trainingData_PCclusters$sw <- cl_silh_width  # add silhouette width to the result

## Save
fname <- paste0(trainingDatasets, "_PCclusters.rds")
saveRDS(trainingData_PCclusters, file.path(out_dir, fname))
