## March 8, 2023 - Bukola Ajanaku
# Creating summarized Heatmap for Habenual progress report!
# qrsh -l mem_free=30G,h_vmem=30G

library("SingleCellExperiment")
library("here")
library("sessioninfo")
library("jaffelab")
library("ggplot2")
library("tidyverse")
library("ComplexHeatmap")
library("spatialLIBD")

load(here("processed-data", "04_snRNA-seq", "sce_objects", 
          "sce_post_09_clustered_qc.Rdata"))

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
    # Astrocyte       Endo Excit.Thal Inhib.Thal      LHb.1      LHb.2      LHb.3 
    # 538         38       1800       7612        201        266        134 
    # LHb.4      LHb.5      LHb.6      LHb.7      MHb.1      MHb.2      MHb.3 
    # 477         83         39       1014        152        540         18 
    # Microglia      Oligo        OPC 
    # 145       2178       1796 

# Pseudobulking to create compressed sce object
## faking the pseudobulking function out by setting sample as all same sample
sce$FakeSample <- "Br1011"
sce$RealSample <- sce$Sample
sce$Sample <- sce$FakeSample

set.seed(20220907) 
sce_simple_pb_snAnno3 <- registration_pseudobulk(sce,  "snAnno3", "Sample")

# grabbing color schemes: grabColors(), max colors 20
source(here("code", "04_snRNA-seq", "sourcing", "color_Scheme_CT.R"))

# list of marker genes 
official_markers = list(
  "Oligo" = c("MOBP"),
  "OPC" = c("PDGFRA"),
  "Micro" = c("CSF1R"),
  "Astro" = c("AQP4"),
  "Endo" = c("ITIH5"),
  'Neuron' = c('SYT1'),
  'Inh_Neuron' = c('GAD1'), 
  'Exc_Neuron' = c('SLC17A6'), 
  "Thal" = c("LYPD6B"),
  'Hb' = c("POU4F1"), 
  "LHb" = c("HTR2C"),
  "MHb" = c("CHAT"),
  "LHb.A " = c("ATP8B1"),
  "LHb.B" = c("CRH"),
  "LHb.C" = c ("ENTHD1"),
  "LHb.D" = c("TLE2"),
  "LHb.E" = c("LINC01619"),
  "LHb.F" = c("TACR3"),
  "LHb.G" = c("AC008619.1"),
  "MHb.A" = c("EXOC1L"), 
  "MHb.B" = c("CHAT"),
  "MHb.C" = c("BHLHE22")
)

namers <- c("Oligo",
"OPC",
"Micro",
"Astro",
"Endo",
'Neuron',
'Inh_Neuron', 
'Exc_Neuron', 
"Thal",
'Hb', 
"LHb",
"MHb",
"LHb.A",
"LHb.B",
"LHb.C",
"LHb.D",
"LHb.E",
"LHb.F",
"LHb.G",
"MHb.A", 
"MHb.B",
"MHb.C")

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
####### PLOTTING ###############################################################
# creating dir for plots
plot_dir <- here("plots", "05_explore_sce", "03_Plot_Investigate_snAnno", 
                 "07_Plots_for_Progress_Report")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

# 

# function to plot heatmap
pseudoHeater <- function(sce, clusterMethod, markerList){
  # sce = pseudobulked sce object
  # clusterMethod = whichever nearest neighbors method you use (as character string)
  
  # Replacing genes with symbols for heatmap (remember, this is pseudobulked data)
  rownames(sce) <- rowData(sce)$Symbol
  
  # renaming rownnames of colData(sce) based on row annotations
  rownames(colData(sce)) <- paste(colData(sce)[, clusterMethod])
  
  # Making data frame of genes we are interested in annd their general classification
  markTable <- as.data.frame(unlist(markerList)) |> 
    rownames_to_column("cellType") |>
    rename(gene = `unlist(markerList)`) |>
    mutate(cellType = gsub("\\d+", "", cellType)) |>
    filter(gene %in% rowData(sce)$Symbol)
  
  # extracting logcounts from sce and subset just genes in markTable
  markerlogcounts <- logcounts(sce[markTable$gene, ])
  
  # getting z scores
  marker_z_score <- scale(t(markerlogcounts))
  # corner(marker_z_score)
  
  # heatmap columns annotation
  column_ha <- HeatmapAnnotation(
    cell_type = markTable$cellType
  )
  
  # grabbing the annotations per cluster from the sce object
  clusterData <- as.data.frame(colData(sce)[,clusterMethod]) 
  names(clusterData) <- "cellType"
  
  # prepping the colors we want
  # for cell type
  num_pal <- length(unique(clusterData$cellType))
  col_pal_ct <- grabColors(num_pal, start = 4)
  names(col_pal_ct) = unique(clusterData$cellType)
  
  # heatmap row annotationn
  row_ha <- rowAnnotation(
    df = clusterData$cellType,
    col = list(cell_type = col_pal_ct)
  )
  
  heatmapped <- Heatmap(marker_z_score,
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        right_annotation = row_ha,
                        top_annotation = column_ha,
                        column_split = factor(markTable$cellType, levels = markTable$cellType) ,
                        column_title_rot = 30,
                        row_order = order(clusterData$cellType, levels = row_namers)
                        )
  
  print(heatmapped)
  
}

# printing 
pdf(here(plot_dir, "Official_Markers_Heatmap_snAnno3_Simple_Pseudobulk.pdf"), width = 12, height = 8)
#  pseudoHeater(sce_simple_pb_snAnno3, "snAnno3", markerList = official_markers)
  print(heatmapped)
dev.off()


# Done.