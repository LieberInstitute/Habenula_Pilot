## March 17, 2023 - Bukola Ajanaku
# Running mean expression and 1 vs All expression methods on bulkTypeSepHb.
# qrsh -l mem_free=20G,h_vmem=20G

library("SingleCellExperiment")
library("here")
library("tidyverse")
library("sessioninfo")
library("DeconvoBuddies")
library("ggplot2")

# loading regular sce object
load(here("processed-data", "99_paper_figs", "sce_objects", 
          "official_final_sce.RDATA"))
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


dim(sce)
# [1] 33848 17031

old_sce <- sce

# dropping noisy OPC 
sce <- sce[, which(sce$OPC_clean == "Yes")]  
  
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

######## PLOTTING #######
## mean expression plots 
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

## hockey stick plots (mean expression vs 1vAll expression)
# adding color group
marker_stats_SepHb$Top25 <- "No"
marker_stats_SepHb[which(marker_stats_SepHb$rank_ratio <= 25),
                   "Top25"] <- "Yes"

pdf(file = here(plot_dir, "marker_stats_SepHb_top5_hockeystick.pdf"))
  
ggplot(marker_stats_SepHb, aes(ratio, std.logFC)) +
  geom_point(size = 0.5, aes(colour = Top25)) +  
  facet_wrap(~cellType.target, scales = "free") 

dev.off()

######## SAVING #######
save(marker_stats_SepHb, mean_ratio_SepHb, findMark_SepHb, file =  here(new_dir, "bulkTypeSepHb_exp.Rdata"))




