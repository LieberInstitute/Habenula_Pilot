## March 10, 2023 - Bukola Ajanaku
# On snAnno2, I will be running mean expression data before and after combining MHb.2 
# and MHb.3. I will not be renaming any clusters other than those that are being changed in order
# to keep organization clear.
# qrsh -l mem_free=20G,h_vmem=20G
# Terminal 4

library("SingleCellExperiment")
library("here")
library("tidyverse")
library("sessioninfo")

# supposed to load from here, but hasn't been saved yet
# load(here("processed-data", "04_snRNA-seq", "sce_objects", "sce_with_snAnno2.RDATA")))

# loading and fixing 
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
          "sce_post_09_clustered_qc.Rdata"))

# creating snAnno2 which will contain the combined MHb groups
sce$snAnno2 <- sce$snAnno

## LHb.6 is actually Endothelial. Total LHb is now 7 from 8.
sce$snAnno2[sce$snAnno2 == "LHb.6"] <- "Endo"

## creating plot_dir for all downstream plots
# for saving plots
plot_dir <- here("plots", "05_explore_sce", "03_Plot_Investigate_snAnno",
                 "05_New_Annotations_meanExpression")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

#### BEFORE COMBINING MHB.2 AND MHB.3 ##########################################
table(sce$snAnno2)


# grabbing mean ration data
mean_ratio_snAnno2 <- get_mean_ratio2(
  sce,
  cellType_col = "snAnno2",
  assay_name = "logcounts",
  add_symbol = TRUE
)

findMark_manAnnno2 <- findMarkers_1vAll(
  sce,
  assay_name = "counts",
  cellType_col = "snAnno2",
  add_symbol = FALSE,
  mod = "~Sample"
)

## For combined manually annotated clusters 
snAnno_marker_stats_new_2 <- left_join(meanRat_snAnno2, findMark_manAnnno2, by = c("gene", "cellType.target"))

pdf(file = here(plot_dir, "Mean_Ratio_Expression_for_new_snAnno2.pdf"))
for (j in levels(as.factor(sce$snAnno2))) {
  message(j)
  
  print(plot_marker_express(sce, 
                             snAnno_marker_stats_new_2, 
                             n_genes = 5,
                             rank_col = "rank_ratio", 
                             anno_col = "anno_ratio",
                             cellType_col = "snAnno2",
                             cell_type = j)
  )
}
dev.off()


#### AFTER COMBINING MHB.2 AND MHB.3 ###########################################
sce$snAnno3 <- sce$snAnno2

# combining MHb.3 with MHb.2
sce$snAnno3[sce$snAnno3 == "MHb.3"] <- "MHb.2"

# sanity check
table(sce$snAnno3)

# 







# saving
save(sce, here("processed-data", "04_snRNA-seq", "sce_objects", "sce_with_snAnno2_and_snAnno3.RDATA")))






# Not done.