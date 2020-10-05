#' This script created MeSH term table used for PCAmodel annotation. The resulting
#' data frame is added in the PCAGenomicSignatures package and you can access it 
#' through:
#' 
#' # dir <- system.file("extdata", package = "PCAGenomicSignatures")
#' # mesh_table <- readRDS(file.path(dir, "MeSH_terms_1399refinebio.rds"))

 


    
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
