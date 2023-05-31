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

# loading sce object
load(here("processed-data", "04_snRNA-seq",  "sce_objects", "sce_final.Rdata"),
     verbose = TRUE)
  # sce_final 

table(sce_final$final_Annotations)
  # Astrocyte       Endo Excit.Thal Inhib.Thal      LHb.1      LHb.2      LHb.3 
  # 538         38       1800       7612        201        266        134 
  # LHb.4      LHb.5      LHb.6      LHb.7      MHb.1      MHb.2      MHb.3 
  # 477         83         39       1014        152        540         18 
  # Microglia      Oligo        OPC 
  # 145       2178       1796 
    # has no Hb cluster 

# creating plot directory
plot_dir <- here("plots", "99_paper_figs", "03_Heatmap_Markers_Summary")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

# sourcing official color palette 
source(file = here("code", "99_paper_figs", "source_colors.R"))
  # bulk_colors and sn_colors

# Pseudobulking to create compressed sce object
## faking the pseudobulking function out by setting sample as all same sample
sce_final$FakeSample <- "Br1011"
sce_final$RealSample <- sce_final$Sample
sce_final$Sample <- sce_final$FakeSample

set.seed(20220907) 
sce_simple_pb_snAnno3 <- registration_pseudobulk(sce_final, "final_Annotations", "Sample")

# list of marker genes 
official_markers = list(
  "Oligo" = c("MOBP"),
  "OPC" = c("PDGFRA"),
  "Micro" = c("CSF1R"),
  "Astro" = c("AQP4"),
  "Endo" = c("ITIH5"),
  "Thal" = c("LYPD6B"),
  "LHb.A" = c("LINC02653"), #  , ATP8B1
  "LHb.B" = c("AC073071.1"),
  "LHb.C" = c ("ENTHD1"),
  "LHb.D" = c("TLE2"),
  "LHb.E" = c("LINC01619"),
  "LHb.F" = c("TACR3"),
  "LHb.G" = c("AC008619.1"),
  "MHb.A" = c("EXOC1L"), 
  "MHb.B" = c("CHAT"),
  "MHb.C" = c("BHLHE22"),
  'Hb' = c("POU4F1"), # BARHL1
  "MHb" = c("CHRNB4"),
  "LHb" = c("HTR2C"),
  'Neuron' = c('SYT1'),
  'Exc_Neuron' = c('SLC17A6'), 
  'Inh_Neuron' = c('GAD1')
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


# explicit color scheme 
color_official_markers = c(
  "Oligo" = c("#4d5802"),
  "OPC"= c("#d3c871"), 
  "Micro" = c("#222222"), 
  "Astro" = c("#8d363c"), 
  "Endo" = c("#ee6c14"), 
  "Thal" = c("#EADDCA"), 
  "LHb.A" = c("#0085af"),
  "LHb.B" = c("#0096FF"),
  "LHb.C" = c ("#89CFF0"), 
  "LHb.D" = c("#6F8FAF"), 
  "LHb.E" = c("#40E0D0"),
  "LHb.F" = c("#008080"), 
  "LHb.G" = c("#7DF9FF"), 
  "MHb.A" = c("#FF00FF"),
  "MHb.B" = c("#FAA0A0"), 
  "MHb.C" = c("#fa246a"),
  "Hb" = c("#702963"), 
  "MHb" = c("#F33A6A"),
  "LHb" = c("#0000FF"), 
  'Neuron' = c('#D8BFD8'),
  'Exc_Neuron' = c("#9e4ad1"), 
  "Inh_Neuron" = c('#b5a2ff')
)


# check colors 
preview_colors <- function(cell_colors) {
  par(las = 2) # make label text perpendicular to axis
  par(mar = c(5, 8, 4, 2)) # increase y-axis margin.
  barplot(rep(1, length(cell_colors)),
          col = cell_colors,
          horiz = TRUE,
          axes = FALSE,
          names.arg = names(cell_colors)
  )
}

png(here(plot_dir, "gene_markers_color_pal.png"))
  preview_colors(color_official_markers)
dev.off()

png(here(plot_dir, "clusters_color_pal.png"))
  preview_colors(sn_colors)
dev.off()

####### PLOTTING ###############################################################
# Plotting ComplexHeatmap
sce = sce_simple_pb_snAnno3
clusterMethod = "final_Annotations"
markerList = official_markers

# Replacing genes with symbols for heatmap (remember, this is pseudobulked data)
rownames(sce) <- rowData(sce)$Symbol

# renaming rownnames of colData(sce) based on row annotations
rownames(colData(sce)) <- paste(colData(sce)[, clusterMethod])

# reordering sce object for plottability
sce_reorder <- sce[unlist(markerList) , row_namers]

# Making data frame of genes we are interested in annd their general classification
markTable <- as.data.frame(unlist(markerList)) |> 
  rownames_to_column("cellType") |>
  rename(gene = `unlist(markerList)`) |>
  mutate(cellType = gsub("\\d+", "", cellType)) |>
  filter(gene %in% rowData(sce_reorder)$Symbol)

# getting z scores
marker_z_score <- scale(t(logcounts(sce_reorder)))
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

# grabbing the annotations per cluster from the sce_reorder object
clusterData <- as.data.frame(colData(sce_reorder)[,clusterMethod]) 
names(clusterData) <- "cellType"

# prepping the colors we want
# for cell type
# num_pal <- length(unique(clusterData$cellType))
# col_pal_ct <- grabColors(num_pal, start = 4)
# names(col_pal_ct) = unique(clusterData$cellType)

# heatmap row annotationn
row_ha <- rowAnnotation(
  Clusters = clusterData$cellType,
  col = list(Clusters = sn_colors)
)

heatmapped <- Heatmap(marker_z_score,
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      right_annotation = row_ha,
                      top_annotation = column_ha,
                      column_split = factor(markTable$cellType, levels = markTable$cellType), 
                      column_title_rot = 30,
                      heatmap_legend_param = list(
                        title = c("Z_Score"),
                        border = "black"
                      ))



# printing 
pdf(here(plot_dir, "Completed_Markers_Heatmap_final_Anno_FINAL.pdf"), width = 12, height = 8)
  heatmapped
dev.off()

#### FOR ONE DRIVE 
# pdf
pdf(here(plot_dir, "forOneDrive", "mfigu_heatmap_progress_report.pdf"), width = 12, height = 8)
  heatmapped
dev.off()

# png
png(here(plot_dir, "forOneDrive", "mfigu_heatmap_progress_report.png"), 
    width = 12, height = 8, units = "in", res = 1200)
heatmapped
dev.off()


# Done.