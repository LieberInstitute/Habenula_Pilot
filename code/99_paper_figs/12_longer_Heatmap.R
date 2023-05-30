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

# making duplicate for later use
sce_final2 <- sce_final
rownames(sce_final2) <- rowData(sce_final2)$Symbol

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
  
  table(sce_final$Sample)
    # Br1011 
    # 17031

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

# grabbing top 5 markers 
marker_genes <- marker_stats |>
  filter(rank_ratio <= 5)

as.list(marker_genes[marker_genes$cellType.target == "LHb", ]$gene)

# list of marker genes 
official_markers = list(
  "Astrocyte" = as.list(marker_genes[marker_genes$cellType.target == "Astrocyte", ]$gene),
  "Endo" = as.list(marker_genes[marker_genes$cellType.target == "Endo", ]$gene),
  "Microglia" = as.list(marker_genes[marker_genes$cellType.target == "Microglia", ]$gene),
  "Oligo" = as.list(marker_genes[marker_genes$cellType.target == "Oligo", ]$gene),
  "OPC" = as.list(marker_genes[marker_genes$cellType.target == "OPC", ]$gene),
  "Thalamus" = as.list(marker_genes[marker_genes$cellType.target == "Thalamus", ]$gene),
  "LHb" = as.list(marker_genes[marker_genes$cellType.target == "LHb", ]$gene),
  "MHb" = as.list(marker_genes[marker_genes$cellType.target == "MHb", ]$gene)
)

set.seed(20220907) 
sce_simple_pb_snAnno3 <- registration_pseudobulk(sce_final2, "final_Annotations", "Sample")

row_namers <- c("Oligo",
                "OPC",
                "Microglia",
                "Astrocyte",
                "Endo",
                "Thalamus",
                "LHb",
                "MHb"
)

# row_namers <- c("Oligo",
#                 "OPC",
#                 "Microglia",
#                 "Astrocyte",
#                 "Endo",
#                 'Inhib.Thal', 
#                 'Excit.Thal', 
#                 "LHb.1",
#                 "LHb.2",
#                 "LHb.3",
#                 "LHb.4",
#                 "LHb.5",
#                 "LHb.6",
#                 "LHb.7",
#                 "MHb.1", 
#                 "MHb.2",
#                 "MHb.3"
# )


# explicit colour scheme for row namers
color_official_markers <- c(
  "Oligo" = c("#4d5802"), 
  "OPC"= c("#d3c871"), 
  "Microglia" = c("#222222"), 
  "Astrocyte" = c("#8d363c"), 
  "Endo" = c("#ee6c14"), 
  "Thalamus"= c("#EADDCA"),
  "LHb" = c("#0085af"),
  "MHb" = c("#fa246a")
)

####### PLOTTING ###############################################################
# Plotting ComplexHeatmap
sce = sce_simple_pb_snAnno3
clusterMethod = "longerAnno"
markerList = official_markers

# Replacing genes with symbols for heatmap (remember, this is pseudobulked data)
rownames(sce) <- rowData(sce)$Symbol

# renaming rownnames of colData(sce) based on row annotations
rownames(colData(sce)) <- paste(colData(sce)[, clusterMethod])

# # reordering sce object for plottability
# sce_reorder <- sce[unlist(markerList) , row_namers]

# Making data frame of genes we are interested in annd their general classification
markTable <- as.data.frame(unlist(markerList)) |> 
  rownames_to_column("cellType") |>
  rename(gene = `unlist(markerList)`) |>
  mutate(cellType = gsub("\\d+", "", cellType)) |>
  filter(gene %in% rowData(sce)$Symbol)

# getting z scores
marker_z_score <- scale(t(logcounts(sce)))
# corner(marker_z_score)

# heatmap columns annotation
column_ha <- HeatmapAnnotation(
  Gene_Marker = markTable$cellType,
  col = list(Gene_Marker = color_official_markers),
  annotation_legend_param = list(
    Gene_Marker = list(
      title = "Gene_Marker" 
    )
  ))


# heatmap row annotationn
row_ha <- rowAnnotation(
  Clusters = sce_final$final_Annotations,
  col = list(Clusters = sn_colors)
)

heatmapped <- Heatmap(marker_z_score,
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      right_annotation = row_ha,
                      top_annotation = column_ha,
                      column_title_rot = 30,
                      heatmap_legend_param = list(
                        title = c("Z_Score"),
                        border = "black"
                      ))



# printing 
pdf(here(plot_dir, "longer_Complex_Heatmap.pdf"), width = 12, height = 8)
  heatmapped
dev.off()

