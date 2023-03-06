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

# for saving plots
plot_dir <- here("plots", "05_explore_sce", "02_Mean_Ratio_Explore")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

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
# in order to run find markers function, Sample must be named donor
sce$donor <- NULL

findMark_wT <- findMarkers_1vAll(
  sce,
  assay_name = "counts",
  cellType_col = "splitSNType",
  add_symbol = FALSE,
  mod = "~Sample"
)


findMark_manAnnno <- findMarkers_1vAll(
  sce,
  assay_name = "counts",
  cellType_col = "snAnno",
  add_symbol = FALSE,
  mod = "~Sample"
)

#### Combining mean ratio with markers #########################################
  ## For uncombined manually annotated clusters
  wT_marker_stats <- left_join(meanRat_wT, findMark_wT, by = c("gene", "cellType.target"))

  ## For combined manually annotated clusters 
  snAnno_marker_stats <- left_join(meanRat_snAnno, findMark_manAnnno, by = c("gene", "cellType.target"))


#### Cell-Type Colors ##########################################################
# gene markers
  
    ## For uncombined manually annotated clusters
#   cell_types <- levels(as.factor(sce$splitSNType))
#   
# pdf(file = here(plot_dir, "cell_type_colors.pdf"))
#   cell_colors <- create_cell_colors(cell_types = cell_types, pallet = "classic", split = NA, preview = TRUE)
# dev.off()
  ## For combined manually annotated clusters 
# cell_types <- levels(as.factor(sce$snAnno))
# 
# pdf(file = here(plot_dir, "cell_type_colors_snAnno.pdf"))
#   cell_colors <- create_cell_colors(cell_types = cell_types, pallet = "classic", split = NA, preview = TRUE)
# dev.off()

#### Plotting Expression #######################################################
# for all 37 clusters
pdf(file = here(plot_dir, "Mean_Ratio_Expression_for_Each_37_Clusters.pdf"))
for (j in levels(as.factor(sce$splitSNType))) {
  message(j)
  
  print(plot_marker_express(sce, 
                    wT_marker_stats, 
                    n_genes = 5,
                    rank_col = "rank_ratio", 
                    anno_col = "anno_ratio",
                    cellType_col = "splitSNType",
                    cell_type = j)
  )
}
dev.off()

# for my grouped clusters
pdf(file = here(plot_dir, "Mean_Ratio_Expression_for_sn_Anno.pdf"))
for (j in levels(as.factor(sce$snAnno))) {
  message(j)
  
  print( plot_marker_express(sce, 
                      snAnno_marker_stats, 
                      n_genes = 5,
                      rank_col = "rank_ratio", 
                      anno_col = "anno_ratio",
                      cellType_col = "snAnno",
                      cell_type = j)
  )
}
dev.off()
  
#### Saving ####################################################################
save(wT_marker_stats, snAnno_marker_stats, findMark_manAnnno, findMark_wT, file =
       here("processed-data", "05_explore_sce", "mean_ratio_data_from_02_mean_ratio_explore_file.Rdata"))
  
  
  
  
  
  
  
  
  
# 
