#' Title: Import RNA-seq data downloaded from refine.bio 
#' 
#' Input: Salmon outputs (`quant.sf` format) of studies with > 20 downloaded samples 
#' - input directory (`in.dir`): refine.bio datasets were downloaded in here
#' - output directory (`out.dir`) : tximport-ed files are saved here
#' Output: tximport output with `.rds` extension
#' Process: Import `_quant.sf` files using `tximport`
#' 
#' Note: Import studies with >20 samples, but for differnet models, different study size 
#' cutoffs were applied later. For example, RAVmodel_536 contains studies with >50 samples.
 



##### For RAVmodel_536 #########################################################
## Working directories
in.dir <- "~/data2/refinebio_download/rna_seq"
# out.dir <- "/nobackup/16tb_b/GenomicSignatures/refinebio/rna_seq"   # 889 studies more than 50 metadata samples/study
out.dir <- "~/data2/refinebio_processed/rna_seq_v2"   # 1579 studies more than 20 downloaded samples/study 

## Load the training dataset information
dir <- system.file("extdata", package = "GenomicSuperSignature")
studyMeta <- read.table(file.path(dir, "studyMeta.tsv.gz"))
ind <- which(studyMeta$downloaded > 20)   # try to import 1579 studies with > 20 downloaded samples
studyNames <- studyMeta$studyName[ind]

## Ensemble annotation database
library(EnsDb.Hsapiens.v86)
txdb <- EnsDb.Hsapiens.v86

##### For all RAVmodels ########################################################
## Inputs required
# wd <- "/mnt/STORE1/16tb_b/refinebio_mice"
# in.dir <- file.path(wd, "rna_seq")
# out.dir <- file.path(wd, "rna_seq_processed")
# studyNames <- c() # study accession numbers as the way it's saved in `in.dir`
# library(EnsDb.Mmusculus.v79) # Change to the correct species
# txdb <- EnsDb.Mmusculus.v79 # Change to the correct species

## Prepare tx2gene file for tximport 
library(tximport)
k <- AnnotationDbi::keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "SYMBOL", "TXNAME")
tx2gene <- tx2gene[,1:2]
tx2gene <- tx2gene[!is.na(tx2gene$SYMBOL),]

## Save combined expression matrix of each study  
for (studyName in studyNames) {
  
  ## collect all the samples
  dir.path <- file.path(in.dir, studyName)   # 6460 studies initially downloaded from refine.bio
  files <- list.files(dir.path)   # all the files in the folder
  files <- files[grepl("_quant.sf", files)]   # only quant.sf files
  files <- file.path(dir.path, files)
  names(files) <- gsub("_quant.sf", "", basename(files))   # assign names
  
  out.path <- file.path(out.dir, studyName)
  if (!dir.exists(out.path)) {dir.create(out.path)}
  
  tryCatch({
      ## tximport
        rseq_counts <- tximport(files = files, type = "salmon", tx2gene = tx2gene,
                                countsFromAbundance = "lengthScaledTPM")
        fname <- paste0(studyName, ".rds")
        saveRDS(rseq_counts, file = file.path(out.path, fname))
  }, error = function(e) {
      # collect studies with error: samples in the same study seem to have different genes
        study_with_error <- data.frame(studyName, conditionMessage(e))
        write.table(study_with_error, 
                    # file = "/nobackup/16tb_b/GenomicSignatures/refinebio/studies_with_error.csv",
                    file = file.path(wd, "studies_with_error.tsv"),
                    append = TRUE, row.names = FALSE, col.names = FALSE)
  })
}