## Global variable for the whole process
main_dir <- "/mnt/STORE1/16tb_b/refinebio_mice"
in_dir <- file.path(main_dir, "rna_seq")
out_dir <- file.path(main_dir, "rna_seq_processed")
wd <- file.path(main_dir, "mouse_colon")
trainingDatasets <- "refinebioRseq_mmu_colon" 
library(dplyr)
library(GenomicSuperSignature)

## 01_Import
studyNames <- c() # study accession numbers as the way it's saved in `in_dir`
library(tximport)
library(EnsDb.Mmusculus.v79)
txdb <- EnsDb.Mmusculus.v79

## 04_Common_Genes
cutoff <- 0.9  # Percent top varying genes to include in RAVmodel

## 05_PCA
n <- 20   # The number PCs to keep

## 06_Clustering
library(factoextra)
d <- 2.25  # for cluster number (line 56)
note <- "29 mouse colon studies"
annotGeneSets <- "MSigDB_mmu_C2"

## 07_Final_Model
library(clusterProfiler)
library(EnrichmentBrowser)
org <- "mmu"
db <- "msigdb"
cat <- "C2"

##### 01_Import ################################################################
## Create directories
if (!dir.exists(out_dir)) {dir.create(out_dir)}
if (!dir.exists(wd)) {dir.create(wd)}

## Prepare tx2gene file for tximport
k <- AnnotationDbi::keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "SYMBOL", "TXNAME")
tx2gene <- tx2gene[,1:2]
tx2gene <- tx2gene[!is.na(tx2gene$SYMBOL),]

## Save combined expression matrix of each study
for (studyName in studyNames) {

    ## collect all the samples
    dir_path <- file.path(in_dir, studyName)
    files <- list.files(dir_path)   # all the files in the folder
    files <- files[grepl("_quant.sf", files)]   # only quant.sf files
    files <- file.path(dir_path, files)
    names(files) <- gsub("_quant.sf", "", basename(files))   # assign names

    out_path <- file.path(out_dir, studyName)
    if (!dir.exists(out_path)) {dir.create(out_path)}

    tryCatch({
        ## tximport
        rseq_counts <- tximport(files = files,
                                type = "salmon",
                                tx2gene = tx2gene,
                                countsFromAbundance = "lengthScaledTPM")
        fname <- paste0(studyName, ".rds")
        saveRDS(rseq_counts, file = file.path(out_path, fname))
    }, error = function(e) {
        ## studies with an error during import
        study_with_error <- data.frame(studyName, conditionMessage(e))
        write.table(study_with_error,
                    file = file.path(wd, "studies_with_error.tsv"),
                    append = TRUE, row.names = FALSE, col.names = FALSE)
    })
}


##### 02_Filtering #############################################################
## Update the available study with successful tximport
import_error <- read.table(file.path(wd, "studies_with_error.tsv"))
studyNames <- studyNames[!studyNames %in% import_error[,1]]

## Log2 transformation
for (studyName in studyNames) {
    dir_path <- file.path(out_dir, studyName)
    file.path <- file.path(dir_path, paste0(studyName, ".rds"))
    countMatrix <- readRDS(file.path)

    # log2(TPM+1) transformation
    count <- countMatrix$counts
    count <- log2(count + 1)

    # Save the results for each studyName in `_count.csv`
    write.table(count, file.path(dir_path, paste0(studyName, "_count.csv")))

    # Standard output to check the progress
    filteringDone <- paste("Successful Filtering done :", studyName)
    print(filteringDone)
}


##### 03_Top_Genes #############################################################
## List containing the ordered list of genes from each study
topGenesInTrainingData <- vector(mode = "list", 
                                 length = length(studyNames)) 

for (studyName in studyNames) {
    dir_path <- file.path(out_dir, studyName)
    x <- data.table::fread(file.path(dir_path, paste0(studyName, "_count.csv")), 
                           showProgress = FALSE)
    
    # order by standard deviation
    y <- apply(x[,-1], 1, sd) 
    names(y) <- x$V1
    y <- y[order(y, decreasing = TRUE)]
    
    topGenesInTrainingData[[studyName]] <- y
    rm(x)
}

print(paste("Finish collecting top", cutoff*100, "% varying genes from", 
            length(studyNames), "studies."))  #--> Standard output

## Save
fname <- paste0("topGenesInTrainingData_", Sys.Date(), ".rds")
saveRDS(topGenesInTrainingData, file.path(wd, fname))


##### 04_Common_Genes ##########################################################
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

topGenes <- topGenes[nzchar(topGenes)] # remove empty character
saveRDS(topGenes, file.path(wd, paste0("topGenes_", length(topGenes),".rds")))


##### 05_PCA ###################################################################
## An empty list for PCA results (rotation and variance)
trainingData_PCA <- vector("list", length(studyNames))
names(trainingData_PCA) <- studyNames

## Calculate sd and mean across all the samples
allSamples <- data.frame(matrix(NA, nrow = length(cg), ncol = 0)) 
rownames(allSamples) <- cg
iter <- 0

for (studyName in studyNames) {
    file_path <- file.path(out_dir, studyName)
    x <- data.table::fread(file.path(file_path, paste0(studyName, "_count.csv")), 
                           showProgress = FALSE)
    x <- data.frame(x[,-1], row.names = x$V1)
    iter <- iter + 1
    allSamples <- cbind(allSamples, x[cg,])
}

print(paste("Calculate standard deviation and mean of", iter, "studies."))
s <- apply(allSamples, 1, sd)
m <- apply(allSamples, 1, mean)
SdMean <- list(sd = s, mean = m)

saveRDS(SdMean, file.path(wd, paste0(trainingDatasets, "_SdMean.rds"))) # sd and mean
saveRDS(allSamples, file.path(wd, paste0(trainingDatasets, ".rds"))) # matrix containing all samples

## Remove non-expressing genes in all samples (m == 0)
non_exp <- which(s == 0) %>% names
s <- s[!names(s) %in% non_exp]
m <- m[!names(m) %in% non_exp]

##### PCA
for (studyName in studyNames) {
    file_path <- file.path(out_dir, studyName)
    x <- data.table::fread(file.path(file_path, paste0(studyName, "_count.csv")), 
                           showProgress = FALSE)
    x <- data.frame(x[,-1], row.names = x$V1)
    x <- x[cg,,drop=FALSE]
    
    # Remove non-expressing genes in all samples (m == 0)
    x <- x[!rownames(x) %in% non_exp,]
    
    # This part will be used if any sample is removed upon filtering
    if (ncol(x) <= 20) {
        final_sample_number <- paste(ncol(x), "samples after filtering.")
        study_with_less_samples <- data.frame(studyName, final_sample_number)
        write.table(study_with_less_samples, 
                    file = file.path(wd, "study_with_less_samples.tsv"),
                    append = TRUE, row.names = FALSE, col.names = FALSE)
        
        print(paste(studyName, "has only", ncol(x), "samples after filtering."))
        next
    }
    
    # Normalization
    x <- sweep(x, 1, m)
    x <- sweep(x, 1, s, "/")
    # x <- preprocessCore::normalize.quantiles(x)   # quantile normalization
    
    # PCA
    pca_res <- prcomp(t(x))   # x is a matrix with genes(row) x samples(column)
    trainingData_PCA[[studyName]]$rotation <- pca_res$rotation[,1:n]
    colnames(trainingData_PCA[[studyName]]$rotation) <- paste0(studyName, ".PC", c(1:n))
    eigs <- pca_res$sdev^2
    pca_summary <- rbind(SD = sqrt(eigs),
                         Variance = eigs/sum(eigs),
                         Cumulative = cumsum(eigs)/sum(eigs))
    trainingData_PCA[[studyName]]$variance <- pca_summary[,1:n]
    colnames(trainingData_PCA[[studyName]]$variance) <- paste0(studyName, ".PC", c(1:n))
    
    rm(x)
}

## Update studyNames after PCA
too_few_samples <- read.table(file.path(wd, "study_with_less_samples.tsv"))
studyNames <- studyNames[!studyNames %in% too_few_samples[,1]]

## Remove empty trainingData_PCA
trainingData_PCA <- trainingData_PCA[names(trainingData_PCA) %in% studyNames]

## Save
fname <- paste0(trainingDatasets, "_PCs_rowNorm.rds")
saveRDS(trainingData_PCA, file.path(wd, fname))


##### 06_Clustering ############################################################
## Combine all PCs
allZ_list <- lapply(trainingData_PCA, function(x) x$rotation)
allZ <- Reduce(cbind, allZ_list)
all <- t(allZ)   # a matrix of PCs (row) x genes (column)
print(paste("Dimension of allZ is", dim(allZ))) 
saveRDS(allZ, file.path(wd, "allZ.rds"))

##### Hierarchical Clustering
## Calculate distance
start <- Sys.time()
res.dist <- factoextra::get_dist(all, method = "spearman")
end <- Sys.time()
t <- end - start
print(paste(t, "took to calculate distance."))  
saveRDS(res.dist, file.path(wd, "res_dist.rds"))

## Cut the tree
k <- round(nrow(all)/d, 0)
start <- Sys.time()
res.hcut <- factoextra::hcut(res.dist, k = k, hc_func = "hclust", 
                             hc_method = "ward.D", hc_metric = "spearman")
end <- Sys.time()
t2 <- end - start
print(paste(round(t2, 2), "seconds took to cut tree."))   
saveRDS(res.hcut, file.path(wd, "res_hcut.rds"))


##### Build avgLoading 
trainingData_PCclusters <- buildAvgLoading(allZ, k, cluster = res.hcut$cluster)

## Silhouette Width
cl <- trainingData_PCclusters$cluster
silh_res <- cluster::silhouette(cl, res.dist)
cl_silh_width <- summary(silh_res)$clus.avg.widths
trainingData_PCclusters$sw <- cl_silh_width  # add silhouette width to the result

## Save
fname <- paste0(trainingDatasets, "_PCclusters.rds")
saveRDS(trainingData_PCclusters, file.path(wd, fname))


##### 07_Final_Model ###########################################################
## Variance Explained from PCA result 
pca_summary <- list()
for (i in seq_along(trainingData_PCA)) {
    pca_summary[[i]] <- trainingData_PCA[[i]]$variance
    names(pca_summary)[i] <- names(trainingData_PCA)[i]
}

##### MeSH terms 
dat_dir <- "~/GSS/GenomicSuperSignaturePaper/inst/extdata" #--> Save in GCP bucket
json_path <- file.path(dat_dir, "sra_to_mesh-000000000000.json") 
x <- jsonlite::stream_in(file(json_path))
mesh_table <- x[x$identifier %in% studyNames,]

len <- length(studyNames)
save(mesh_table, 
     file = file.path(wd, paste0("MeSH_terms_", len, "refinebio.rda")))

## Create bagOfWords and MeSH_freq
MeSH_freq <- table(mesh_table$name) %>% 
    sort(., decreasing = TRUE)  # frequency of each term
for (i in seq_len(nrow(mesh_table))) { 
    mesh_table$bagOfWords[i] <- MeSH_freq[mesh_table$name[i]]
} # add frequency to the main table 

unique_id <- unique(mesh_table$identifier)

## Split MeSH term table to a list of each study, `all_MeSH`
all_MeSH <- vector("list", length(unique_id))
names(all_MeSH) <- unique_id
for (study in unique_id) {
    ind <- grepl(study, mesh_table$identifier)
    mesh_meta <- c("score", "identifier", "concept", "name",
                   "explanation", "bagOfWords") # Currently excludes 'concept'
    all_MeSH[[study]] <- mesh_table[ind, mesh_meta]
}

## In case MeSH term information is not available for all the studies
trainingData_MeSH <- all_MeSH[studyNames] # Subset only to the studies used in the model

##### Build PCAGenomicSignatures object
## Assemble training data
trainDat <- as.data.frame(matrix(nrow = length(pca_summary), ncol = 0))
trainDat$PCAsummary <- pca_summary
row.names(trainDat) <- names(pca_summary)

RAVindex <- as.matrix(trainingData_PCclusters$avgLoading)
colnames(RAVindex) <- paste0("RAV", seq_len(ncol(RAVindex)))

## Construct RAVmodel
RAVmodel <- PCAGenomicSignatures(assays = list(RAVindex = RAVindex), 
                                 trainingData = DataFrame(trainDat))

## metadata
metadata(RAVmodel) <- trainingData_PCclusters[c("cluster", "size", "k", "n")]
names(metadata(RAVmodel)$size) <- paste0("RAV", seq_len(ncol(RAVmodel)))
geneSets(RAVmodel) <- annotGeneSets
metadata(RAVmodel)$MeSH_freq <- MeSH_freq
updateNote(RAVmodel) <- note
metadata(RAVmodel)$version <- "0.1.0"

## colData
colData(RAVmodel)$RAV <- paste0("RAV", seq_len(ncol(RAVmodel)))
studies(RAVmodel) <- trainingData_PCclusters$studies
silhouetteWidth(RAVmodel) <- trainingData_PCclusters$sw

## trainingData
mesh(RAVmodel) <- trainingData_MeSH[rownames(trainingData(RAVmodel))]

## Save pre-GSEA version
fname <- paste0(trainingDatasets, "_RAVmodel_", 
                format(Sys.Date(), format="%Y%m%d"), ".rds")
saveRDS(RAVmodel, file.path(wd, fname))


##### GSEA #####################################################################
gsea_dir <- file.path(wd, paste0("gsea_", annotGeneSets))
if(!dir.exists(gsea_dir)) {dir.create(gsea_dir)}

gs <- EnrichmentBrowser::getGenesets(org = org, 
                                     db = db, 
                                     cat = cat, # category 
                                     gene.id.type = "SYMBOL")
term2gene <- stack(gs)
colnames(term2gene) <- c("entrez_gene", "gs_name")
term2gene <- term2gene[,c(2,1)] # the order of columns should be term and gene.

## From `model_building/inst/scripts/build_gsea_DB.R` script
for (i in seq_len(ncol(RAVmodel))) {
    fname <- paste0("gsea_", i, ".rds")
    fpath <- file.path(gsea_dir, fname)
    
    geneList <- RAVindex(RAVmodel)[,i]
    geneList <- sort(geneList, decreasing = TRUE)
    res <- clusterProfiler::GSEA(geneList, TERM2GENE = term2gene,
                                 pvalueCutoff = 0.05, seed = TRUE)
    saveRDS(res, fpath)
}

## Load searchPathways function
func_path <- system.file("scripts/searchPathways.R", package = "GenomicSuperSignature")
source(func_path)

gsea_all <- searchPathways(RAVmodel, gsea_dir)  
gsea(RAVmodel) <- gsea_all

## Save post-GSEA version
fname <- paste0(trainingDatasets, "_RAVmodel_", annotGeneSets, 
                "_", format(Sys.Date(), format="%Y%m%d"),".rds")
saveRDS(RAVmodel, file.path(wd, fname))
