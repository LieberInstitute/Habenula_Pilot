## March 17, 2023 - Bukola Ajanaku
# Running mean expression and 1 vs All expression methods on bulkTypeSepHb.
# qrsh -l mem_free=20G,h_vmem=20G

library("SingleCellExperiment")
library("here")
library("tidyverse")
library("sessioninfo")
library("DeconvoBuddies")

# loading regular sce object
load(here("processed-data", "06_deconvolution", "sce_objects", "sce_first_bulkTypes.RDATA"))

# creating directories for saving  
  ## output save directory
  new_dir <- here("processed-data", "06_deconvolution", "sce_objects", "expression_data")
  if(!dir.exists(new_dir)){
    dir.create(new_dir)
  }
  
  ## plotting directory
  plot_dir <- here("plots", "06_deconvolution", "expression_data")
  if(!dir.exists(plot_dir)){
    dir.create(plot_dir)
  }

##### Running expression code! #####
    # mean ratio
  mean_ratio_SepHb <- get_mean_ratio2(
    sce,
    cellType_col = "bulkTypeSepHb",
    assay_name = "logcounts",
    add_symbol = TRUE
  )

    # 1vAll Expression
  findMark_SepHb <- findMarkers_1vAll(
    sce,
    assay_name = "counts",
    cellType_col = "bulkTypeSepHb",
    add_symbol = FALSE,
    mod = "~Sample"
  )

## For combined manually annotated clusters 
marker_stats_SepHb <- left_join(mean_ratio_SepHb, findMark_SepHb, by = c("gene", "cellType.target"), 
                                multiple = "all")

pdf(file = here(plot_dir, "marker_stats_SepHb_top5_meanExp.pdf"))
for (j in levels(as.factor(sce$bulkTypeSepHb))) {
  message(j)
  
  print(plot_marker_express(sce, 
                            marker_stats_SepHb, 
                            n_genes = 5,
                            rank_col = "rank_ratio", 
                            anno_col = "anno_ratio",
                            cellType_col = "bulkTypeSepHb",
                            cell_type = j)
  )
}
dev.off()

# Saving 
save(marker_stats_SepHb, mean_ratio_SepHb, findMark_SepHb, file =  here(new_dir, "bulkTypeSepHb_exp.Rdata"))




