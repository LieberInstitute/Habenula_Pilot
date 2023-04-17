## March 17, 2023 - Bukola Ajanaku
# Running mean expression and 1 vs All expression methods on bulkTypeBroadHb.
# qrsh -l mem_free=20G,h_vmem=20G

library("SingleCellExperiment")
library("here")
library("tidyverse")
library("sessioninfo")
library("DeconvoBuddies")
library("ggplot2")

# loading regular sce object
load(here("processed-data", "06_deconvolution", "sce_objects", "sce_first_bulkTypes.RDATA"),
     verbose = TRUE)

# creating directories for saving  
## output save directory
new_dir <- here("processed-data", "06_deconvolution", "sce_objects", "expression_data")
## plotting directory
plot_dir <- here("plots", "06_deconvolution", "expression_data")


##### Running expression code! #####
# mean ratio
mean_ratio_BroadHb <- get_mean_ratio2(
  sce,
  cellType_col = "bulkTypeBroadHb",
  assay_name = "logcounts",
  add_symbol = TRUE
)

# 1vAll Expression
findMark_BroadHb <- findMarkers_1vAll(
  sce,
  assay_name = "counts",
  cellType_col = "bulkTypeBroadHb",
  add_symbol = FALSE,
  mod = "~Sample"
)

## For combined manually annotated clusters 
marker_stats_BroadHb <- left_join(mean_ratio_BroadHb, findMark_BroadHb, by = c("gene", "cellType.target"), 
                                multiple = "all")

######## PLOTTING #######
## mean expression plots 
pdf(file = here(plot_dir, "marker_stats_BroadHb_top5_meanExp.pdf"))
for (j in levels(as.factor(sce$bulkTypeBroadHb))) {
  message(j)
  
  print(plot_marker_express(sce, 
                            marker_stats_BroadHb, 
                            n_genes = 5,
                            rank_col = "rank_ratio", 
                            anno_col = "anno_ratio",
                            cellType_col = "bulkTypeBroadHb",
                            cell_type = j)
  )
}
dev.off()

## hockey stick plots (mean expression vs 1vAll expression)
# adding color group
marker_stats_BroadHb$Top25 <- "No"
marker_stats_BroadHb[which(marker_stats_BroadHb$rank_ratio <= 25),
                   "Top25"] <- "Yes"

pdf(file = here(plot_dir, "marker_stats_BroadHb_top5_hockeystick.pdf"))

ggplot(marker_stats_BroadHb, aes(ratio, std.logFC)) +
  geom_point(size = 0.5, aes(colour = Top25)) +  
  facet_wrap(~cellType.target, scales = "free") 

dev.off()

######## SAVING #######
save(marker_stats_BroadHb, mean_ratio_BroadHb, findMark_BroadHb, file =  here(new_dir, "bulkTypeBroadHb_exp.Rdata"))




