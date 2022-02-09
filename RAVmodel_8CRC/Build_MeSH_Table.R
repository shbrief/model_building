#' This script created MeSH term table used for RAVmodel_CRC annotation. 

## Load the training dataset information
# dir <- system.file("extdata", package = "GenomicSuperSignature")
# studyMeta <- read.table(file.path(dir, "studyMeta.tsv.gz"))
# ind <- which(studyMeta$imported == TRUE)
# allStudies <- studyMeta$studyName[ind]

allStudies <- c("GSE13294")

##### MeSH #####################################################################
dat_dir <- "~/data2/GenomicSuperSignaturePaper/inst/extdata"
json_path <- file.path(dat_dir, "sra_to_mesh-000000000000.json")
x <- jsonlite::stream_in(file(json_path))
x <- x[x$identifier %in% allStudies,]
mesh_table <- x
save(mesh_table, file = file.path(dat_dir, "MeSH_terms_1399refinebio.rda"))
