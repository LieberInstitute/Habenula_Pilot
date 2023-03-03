## March 3, 2023 - Bukola Ajanaku
# Running mean ratio expression analysis to help confirm annotated identity of 
# kept clusters (all 37).
# Inspired by: http://research.libd.org/DeconvoBuddies/index.html
# qrsh -l mem_free=20G,h_vmem=20G

library("SingleCellExperiment")
library("here")
library("sessioninfo")
library("DeconvoBuddies")
# library("jaffelab")
library("tidyverse")

# loading sce with annotation names 
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
          "sce_post_09_clustered_qc.Rdata"))

#### Mean Ratios ###############################################################
meanRat_wT <- get_mean_ratio2(
  sce,
  cellType_col = "splitSNType",
  assay_name = "logcounts",
  add_symbol = TRUE
)

meanRat_snAnno <- get_mean_ratio2(
  sce,
  cellType_col = "snAnno",
  assay_name = "logcounts",
  add_symbol = TRUE
)

#### Find Markers ##############################################################
sce$donor <- sce$Sample 

findMark_wT <- findMarkers_1vAll(
  sce,
  assay_name = "counts",
  cellType_col = "splitSNType",
  add_symbol = FALSE
)


findMark_manAnnno <- findMarkers_1vAll(
  sce,
  assay_name = "counts",
  cellType_col = "snAnno",
  add_symbol = FALSE
)

#### Combining mean ratio with markers #########################################
  ## For uncombined manually annotated clusters
  wT_marker_stats <- left_join(meanRat_wT, findMark_wT, by = c("gene", "cellType.target"))

  ## For combined manually annotated clusters 
  snAnno_marker_stats <- left_join(meanRat_snAnno, findMark_manAnnno, by = c("gene", "cellType.target"))


#### Cell-Type Colors ##########################################################
  


  
  
#### Saving ####################################################################
save(wT_marker_stats, snAnno_marker_stats, findMark_manAnnno, findMark_wT, file =
       here("processed-data", "05_explore_sce", "mean_ratio_data_from_02_mean_ratio_explore_file.Rdata"))
  
  
  
  
  
  
  
  
  
# 
