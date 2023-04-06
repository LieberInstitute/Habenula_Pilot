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
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
          "sce_post_09_clustered_qc.Rdata"))

# creating dir for plots
plot_dir <- here("plots", "99_paper_figs", "03_Heatmap_Markers_Summary")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

# creating snAnno3
sce$snAnno3 <- sce$snAnno

# renaming LHb.6 to Endo
sce$snAnno3[sce$snAnno3 == "LHb.6"] <- "Endo"
# combining MHb.3 with MHb.2
sce$snAnno3[sce$snAnno3 == "MHb.3"] <- "MHb.2"
# changing names for specific clusters
sce$snAnno3[sce$snAnno3 == "LHb.7"] <- "LHb.6"
sce$snAnno3[sce$snAnno3 == "LHb.8"] <- "LHb.7"
sce$snAnno3[sce$snAnno3 == "MHb.4"] <- "MHb.3"

# dropping confusing Hb cluster
sce <- sce[ , which(sce$snAnno3 != "Hb")]

table(sce$snAnno3)
# Astrocyte       Endo Excit.Thal Inhib.Thal  LHb.1      LHb.2      LHb.3 
# 538         38       1800       7612        201        266        134 
# LHb.4      LHb.5      LHb.6      LHb.7      MHb.1      MHb.2      MHb.3 
# 477         83         39       1014        152        540         18 
# Microglia   Oligo      OPC 
# 145         2178       1796 

# Pseudobulking to create compressed sce object
## faking the pseudobulking function out by setting sample as all same sample
sce$FakeSample <- "Br1011"
sce$RealSample <- sce$Sample
sce$Sample <- sce$FakeSample

set.seed(20220907) 
sce_simple_pb_snAnno3 <- registration_pseudobulk(sce,  "snAnno3", "Sample")


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
## # glia cells: greens (5), thalamus (yellow #d6c68c), habenula (lateral  vs medial), 
# neurons: purple (inhib (red) vs excit (blue))

color_official_markers = c(
  "Oligo" = c("#5C4033"), # dark grown
  "OPC"= c("#899499"), # pewter (grey)
  "Micro" = c("#4CBB17"), # kelly green
  "Astro" = c("#CC5500"), # burnt orange
  "Endo" = c("#702963"), # byzantium
  "Thal" = c("#FF69B4"), # hot pink
  "LHb.A" = c("#5F9EA0"), # cadet Blue
  "LHb.B" = c("#5D3FD3"), # iris
  "LHb.C" = c ("#4682B4"), #Steel Blue
  "LHb.D" = c("#1F51FF"), # neon blue
  "LHb.E" = c("#6495ED"), # Cornflower Blue
  "LHb.F" = c("#088F8F"), # Blue Green 
  "LHb.G" = c("#4169E1"), # royal bluee
  "MHb.A" = c("#40E0D0"), # Turquoise
  "MHb.B" = c("#96DED1"), #Robin Egg Blue
  "MHb.C" = c("#7DF9FF"), # Electric Blue
  "Hb" = c("#6082B6"), # glaucous (ashy blue like denim)
  "MHb" = c("#00FFFF"), # aqua (bright blue)
  "LHb" = c("#3F00FF"), # indigo (deeper blue)
  'Neuron' = c('#800080'), # purple
  'Exc_Neuron' = c('#0818A8'), # zaffre (dark blue)
  "Inh_Neuron" = c('#FF0000') #  red
)
#833866

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


## glia cells: glia (5 colors), thalamus (2 reds), 
color_row_namers <- c( "Oligo" = c("#5C4033"), # dark grown
                       "OPC"= c("#899499"), # pewter (grey)
                       "Microglia" = c("#4CBB17"), # kelly green
                       "Astrocyte" = c("#CC5500"), # burnt orange
                       "Endo" = c("#702963"), # byzantium
                       "Inhib.Thal" = c('#FF0000'), #red
                       "Excit.Thal" = c('#0818A8'), # zaffre (dark blue)
                       "LHb.1" = c("#5F9EA0"), # cadet Blue
                       "LHb.2" = c("#5D3FD3"), # iris
                       "LHb.3" = c ("#4682B4"), #Steel Blue
                       "LHb.4" = c("#1F51FF"), # neon blue
                       "LHb.5" = c("#6495ED"), # Cornflower Blue
                       "LHb.6" = c("#088F8F"), # Blue Green 
                       "LHb.7" = c("#4169E1"), # royal bluee
                       "MHb.1" = c("#40E0D0"), # Turquoise
                       "MHb.2" = c("#96DED1"), #Robin Egg Blue
                       "MHb.3" = c("#7DF9FF") # light blue
)

png(here(plot_dir, "clusters_color_pal.png"))
  preview_colors(color_row_namers)
dev.off()

####### PLOTTING ###############################################################
# Plotting ComplexHeatmap
sce = sce_simple_pb_snAnno3
clusterMethod = "snAnno3"
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
  col = list(Clusters = color_row_namers)
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
pdf(here(plot_dir, "Completed_Markers_Heatmap_snAnno3_Simple_Pseudobulk.pdf"), width = 12, height = 8)
  heatmapped
dev.off()


# Done.