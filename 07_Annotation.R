#' Title: Annotating PCclusters using MeSH and GSEA
#' Input: 
#' Output: 
#' 
#' Process:





## Working directories
out.dir <- "~/data2/PCAGenomicSignatureLibrary/refinebioRseq/PCAmodel_536"
    
## Load the training dataset information
dir <- system.file("extdata", package = "PCAGenomicSignatures")
studyMeta <- read.table(file.path(dir, "studyMeta.tsv"))
ind <- which(studyMeta$PCAmodel_536 == TRUE)
allStudies <- studyMeta$studyName[ind]

##### MeSH #####################################################################
json_path <- "~/data2/Genomic_Super_Signature/data/sra_to_mesh-000000000000.json"
x <- jsonlite::stream_in(file(json_path))
x <- x[x$identifier %in% allStudies,]
saveRDS(x, "~/data2/PCAGenomicSignatures/inst/extdata/MeSH_terms_1399refinebio.rds")

##### GSEA #####################################################################
dir2 <- system.file("script", package = "PCAGenomicSignatures")
gsea_script <- file.path(dir2, "build_gsea_DB.R")
source(gsea_script)  # This is the processing script. Doule-check the details.

searchPathways_func <- file.path(dir2, "searchPathways.R")
source(searchPathways_func)  # load the function

gsea_all <- searchPathways(PCAmodel, 
                           file.path(out.dir, "gsea"))  # GSEA DB is saved here
gsea(PCAmodel) <- gsea_all
updateNote(PCAmodel) <- "536 refine.bio studies/ top 90% varying genes/ GSEA with MSigDB C2.all"
