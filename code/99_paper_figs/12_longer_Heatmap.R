## April 6, 2023 - Bukola Ajanaku
# This is the exact code I used for the original progress report. Copy and pasted!
# Simply changed plotting directory!!!
# qrsh -l mem_free=30G,h_vmem=30G

library("SingleCellExperiment")
library("here")
library("sessioninfo")
library("jaffelab")
library("ggplot2")
library("tidyverse")
library("ComplexHeatmap")
library("spatialLIBD")
library("DeconvoBuddies")

# loading sce object
load(here("processed-data", "04_snRNA-seq",  "sce_objects", "sce_final.Rdata"),
     verbose = TRUE)
# sce_final 

# creating plot directory
plot_dir <- here("plots", "99_paper_figs", "12_longer_Heatmap")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

# sourcing official color palette 
source(file = here("code", "99_paper_figs", "source_colors.R"))
# bulk_colors and sn_colors

# messing with columns for registration pseudobulk
table(sce_final$RealSample)
  # Br1092 Br1204 Br1469 Br1735 Br5555 Br5558 Br5639 
  # 3208   1193   2150   3101   3542    778   3059

table(sce_final$final_Annotations)
  # grabbing the resolution we need
sce_final$longerAnno <- sce_final$final_Annotations
# making separated Hb (2)
sce_final$longerAnno[sce_final$longerAnno %in% grep("^LHb\\.", unique(sce_final$longerAnno), value = TRUE)] <-  'LHb'
sce_final$longerAnno[sce_final$longerAnno %in% grep("^MHb\\.", unique(sce_final$longerAnno), value = TRUE)] <-  'MHb'
sce_final$longerAnno[sce_final$longerAnno %in% grep("Thal*", unique(sce_final$longerAnno), value = TRUE)] <- "Thalamus"

table(sce_final$longerAnno)
    # Astrocyte      Endo       LHb       MHb Microglia     Oligo       OPC  Thalamus 
    # 538        38      2214       710       145      2178      1796      9412 

set.seed(20220907) 
sce_simple_pb_snAnno3 <- registration_pseudobulk(sce_final, "longerAnno", "Sample")

# finding necessary gene markers 
# adding symbols
rownames(sce_final) <- rowData(sce_final)$Symbol

# pre-bisque
# Creating mean_ratios based on our specified annotations
ratios <- get_mean_ratio2(sce_final,
                          cellType_col = "longerAnno",
                          assay_name = "logcounts",
                          add_symbol = TRUE)

# Using the 1 vs All standard fold change for each gene x cell type
fc <- findMarkers_1vAll(sce_final,
                        assay_name = "counts",
                        cellType_col = "longerAnno",
                        add_symbol = FALSE,
                        mod = "~RealSample",
                        verbose = TRUE
)

# combining data
marker_stats <- left_join(ratios, fc, by = c("gene", "cellType.target"))



# list of marker genes 
official_markers = list(
  "Astrocyte" = c("MOBP"),
  "Endo" = c(),
  "Microglia" = c(),
  "Oligo" = c(),
  "OPC" = c(),
  "Thalamus" = c(),
  "LHb" = c(),
  "MHb" = c()
)


row_namers <- c("Oligo",
                "OPC",
                "Microglia",
                "Astrocyte",
                "Endo",
                'Inhib.Thal', 
                'Excit.Thal', 
                "LHb.1",
                "LHb.2",
                "LHb.3",
                "LHb.4",
                "LHb.5",
                "LHb.6",
                "LHb.7",
                "MHb.1", 
                "MHb.2",
                "MHb.3"
)
