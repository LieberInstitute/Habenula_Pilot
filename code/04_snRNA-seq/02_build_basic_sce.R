# January 6, 2023 - Bukola Ajanaku
# Building basic sce objects. This follows the initial pca analysis for bulk 
# Habenula data. Only 7 samples, all control.
# Based on: 
# 1) https://github.com/LieberInstitute/Habenula_Bulk/blob/master/code/08_snRNA-seq_Erik/20210323_human_hb_neun.R
# 2) https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/code/03_build_sce/01_build_basic_sce.R
# 3) https://github.com/LieberInstitute/Habenula_Bulk/blob/master/code/09_snRNA-seq_re-processed/00_get_droplet_scores.R
# qrsh -l mem_free=75G,h_vmem=75G

# Loading relevant libraries
library("SingleCellExperiment")
library("SpatialExperiment")
library("DropletUtils")
library("here")
library("rtracklayer")
library("lobstr")
library("sessioninfo")
library("dplyr")
library("scater")
library("purrr")

# Sample names 
sample_list <- list.files(here("processed-data", "07_cellranger"))

sce_list <- map(sample_list, function(sample){
  fn10x <-  here("processed-data", "07_cellranger", sample, "outs", "raw_feature_bc_matrix") 
  fnDropScore <- here("processed-data", "04_snRNA-seq", "01_get_droplet_scores", paste0("droplet_scores_", sample, ".Rdata"))
  
  sce <- read10xCounts(fn10x, col.names=TRUE)
  load(fnDropScore, verbose = TRUE)
  ncol_preDrop = ncol(sce)
  # 1148322
  
  sce <- sce[, which(e.out$FDR <= 0.001)]
  ncol_postDrop = ncol(sce)
  # 3622
  
  message(sample, ": Pre-Drop = ", ncol_preDrop, " and Post-Drop = ", ncol_postDrop)

  return(sce)
})

# Match rownames and appending by column
make_Rows <- rownames(sce_list[[1]])
sce_list <- map(sce_list, ~.x[make_Rows,])
# identical(rownames(sce_list2[[1]]), rownames(sce_list2[[7]]))
# TRUE

# OFFICIAL SCE HB OBJECT PRE-QC. 
sce_hb_preQC <- do.call("cbind", sce_list)

save(sce_hb_preQC, file = here("processed-data", "04_snRNA-seq", "sce_objects", 
                               "sce_hb_preQC.Rdata"))
# sgejobs::job_single('02_build_basic_sce', create_shell = TRUE, queue= 'bluejay', memory = '100G', command = "Rscript 02_build_basic_sce.R")


## QUALITY CONTROL
# Thought process: (this is 10x Genomics data)
# 1) Based on OSCA book, can use mean absolute deviation (MADs) approach to 
# determine outliers but this is not as "straightforward" as it may through 
# away neurons.
# 2) There is a test using UMI/barcode rank (the knee plots), to determine 
# individual thresholds for what is too low of quality but this may drop cells 
# with organically low RNA content.
# 3) This is a single-nuc RNA project meaning that dropping high mito content
# is a useful QC metric because it shows samples where cytoplasm 
# was not fully or successfully stripped.
# Game plan:
# 1) Drop empty droplets. (script 1, completed)
# 2) Build sce objects. (this script, script 2)
# 3) Dropping high mito content (script 3, next script)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
