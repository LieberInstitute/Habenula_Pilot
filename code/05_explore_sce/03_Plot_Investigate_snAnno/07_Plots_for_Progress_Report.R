## March 8, 2023 - Bukola Ajanaku
# Creating summarized Heatmap for Habenual progress report!
# qrsh -l mem_free=20G,h_vmem=20G

library("SingleCellExperiment")
library("here")
library("sessioninfo")
library("jaffelab")
library("ggplot2")
library("tidyverse")
library("ComplexHeatmap")

# loading pseudobulked data (snAnno3 with combined MHb.2 and MHb.3)
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
          "sce_new_pseudobulk_with_snAnno_version2and3.Rdata"), verbose = TRUE)
sce <- sce_new_pb_snAnno3

# grabbing color schemes: grabColors(), max colors 20
source(here("code", "04_snRNA-seq", "sourcing", "color_Scheme_CT.R"))

# changing names for specific clusters
sce$snAnno3[sce$snAnno3 == "Hb"] <- "Neuron.Ambig"
sce$snAnno3[sce$snAnno3 == "LHb.7"] <- "LHb.6"
sce$snAnno3[sce$snAnno3 == "LHb.8"] <- "LHb.7"
sce$snAnno3[sce$snAnno3 == "MHb.4"] <- "LHb.3"

table(sce$snAnno3)
  # Astrocyte         Endo   Excit.Thal   Inhib.Thal        LHb.1        LHb.2 
  # 3            2            6            6            4            4 
  # LHb.3        LHb.4        LHb.5        LHb.6        LHb.7        MHb.1 
  # 5            4            3            2            5            3 
  # MHb.2    Microglia Neuron.Ambig        Oligo          OPC 
  # 4            3            1            4            6 

# list of marker genes 
progress.markers = list(
  'Neuron' = c('SYT1', 'SNAP25', "SYT4", "SYP"), 
  'Ihib Neuron' = c('GAD1', 'GAD2'), 
  'Excit Neuron' = c('SLC17A6', 'SLC17A7'), 
  'Classic Habenula' = c("POU4F1",'CHRNB3', 'LMO3', "LRRC55", "MAPK4", "SSTR2", "IRX2", 
                             "NR4A2", "VAV2", "ADCYAP1", "RPRM", "CHRNA3"), 
  "Classic LHb Neuron" = c("HTR2C", "MMRN1", "GAP43", "PARM1", "PLCH1"),
  "Classic MHb Neuron " = c("NEUROD1", "SLC17A7", "TAC3", "CHAT"),
  "MHb.N Group 1" = c("SPON1", "SATB1", "FABP5"),
  "MHb.N Group 2" = c("CCK", "TAC1"),
  "LHb.N Group 1" = c("CHRM3", "VGF", "GABRA1"),
  "LHb.N Group 2" = c("GPR151"),
  "Thal.Med" = c("ADARB2", "LYPD6", "LYPD6B", "EPHA4", "SULF1", "LINC02253"),
  "Endo/Mural" = c("CARMN", "NOTCH3", "SLC38A11"),
  "Oligo" = c("MOBP", "MBP"),
  "OPC" = c("PDGFRA", "VCAN"),
  "Microglia" = c("C3", "CSF1R"),
  "Astrocyte" = c("GFAP", "AQP4", "ETNPPL"),
  "Main LHb Cluster Markers" = c("LINC02653", "SYT2", "AC073071.1", "CRH", "ENTHD1",
                        "GRAP2", "TLE2", "CBLN1", "LINC01619", "CCND2",
                        "TACR3", "ESRP1", "AC008619.1", "LINC02296"),
  "Main MHb Cluster Markers" = c("EXOC1L", "F13A1", "CHAT", "AC079760.2", "BHLHE22",
                                 "PKP2")
)

# list of marker genes 
progress.markers_split = list(
  'Neuron' = c('SYT1', 'SNAP25', "SYT4", "SYP"), 
  'Ihib Neuron' = c('GAD1', 'GAD2'), 
  'Excit Neuron' = c('SLC17A6', 'SLC17A7'), 
  'Classic Habenula' = c("POU4F1",'CHRNB3', 'LMO3', "LRRC55", "MAPK4", "SSTR2", "IRX2", 
                         "NR4A2", "VAV2", "ADCYAP1", "RPRM", "CHRNA3"), 
  "Classic LHb Neuron" = c("HTR2C", "MMRN1", "GAP43", "PARM1", "PLCH1"),
  "Classic MHb Neuron " = c("NEUROD1", "SLC17A7", "TAC3", "CHAT"),
  "MHb.N Group 1" = c("SPON1", "SATB1", "FABP5"),
  "MHb.N Group 2" = c("CCK", "TAC1"),
  "LHb.N Group 1" = c("CHRM3", "VGF", "GABRA1"),
  "LHb.N Group 2" = c("GPR151"),
  "Thal.Med" = c("ADARB2", "LYPD6", "LYPD6B", "EPHA4", "SULF1", "LINC02253"),
  "Endo/Mural" = c("CARMN", "NOTCH3", "SLC38A11"),
  "Oligo" = c("MOBP", "MBP"),
  "OPC" = c("PDGFRA", "VCAN"),
  "Microglia" = c("C3", "CSF1R"),
  "Astrocyte" = c("GFAP", "AQP4", "ETNPPL"),
  "LHb C " = c("LINC02653", "SYT2"),
  "LHb Cl" = c("AC073071.1", "CRH"),
  "LHb Clu" = c ("ENTHD1","GRAP2"),
  "LHb Clus" = c("TLE2", "CBLN1"),
  "LHb Clust" = c("LINC01619", "CCND2"),
  "LHb Cluste" = c("TACR3", "ESRP1"),
  "LHb Cluster" = c("AC008619.1", "LINC02296"),
  "MHb Cl" = c("EXOC1L", "F13A1"), 
  "MHb Clus" = c("CHAT", "AC079760.2"),
  "MHb Cluste" = c("BHLHE22", "PKP2")
)

# list of marker genes 
single.markers_split = list(
  'Neuron' = c('SYT1'), 
  'Ihib Neuron' = c('GAD1'), 
  'Excit Neuron' = c('SLC17A6'), 
  'Classic Habenula' = c("POU4F1"), 
  "Classic LHb Neuron" = c("HTR2C"),
  "Classic MHb Neuron " = c("NEUROD1"),
  "MHb.N Group 1" = c("SPON1"),
  "MHb.N Group 2" = c("CCK"),
  "LHb.N Group 1" = c("CHRM3"),
  "LHb.N Group 2" = c("GPR151"),
  "Thal.Med" = c("ADARB2"),
  "Oligo" = c("MOBP"),
  "OPC" = c("PDGFRA"),
  "Microglia" = c("C3"),
  "Astrocyte" = c("GFAP"),
  "LHb C " = c("LINC02653"),
  "LHb Cl" = c("AC073071.1"),
  "LHb Clu" = c ("ENTHD1"),
  "LHb Clus" = c("TLE2"),
  "LHb Clust" = c("LINC01619"),
  "LHb Cluste" = c("TACR3"),
  "LHb Cluster" = c("AC008619.1"),
  "MHb Cl" = c("EXOC1L"), 
  "MHb Clus" = c("CHAT"),
  "MHb Cluste" = c("BHLHE22")
)

# list of marker genes 
single.markers_split_simplified = list(
  'Classic Habenula' = c("POU4F1"), 
  "Classic LHb Neuron" = c("HTR2C"),
  "Classic MHb Neuron " = c("NEUROD1"),
  "Thal.Med" = c("ADARB2"),
  "Oligo" = c("MOBP"),
  "OPC" = c("PDGFRA"),
  "Microglia" = c("C3"),
  "Astrocyte" = c("GFAP"),
  "LHb C " = c("LINC02653"),
  "LHb Cl" = c("AC073071.1"),
  "LHb Clu" = c ("ENTHD1"),
  "LHb Clus" = c("TLE2"),
  "LHb Clust" = c("LINC01619"),
  "LHb Cluste" = c("TACR3"),
  "LHb Cluster" = c("AC008619.1"),
  "MHb Cl" = c("EXOC1L"), 
  "MHb Clus" = c("CHAT"),
  "MHb Cluste" = c("BHLHE22")
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
  clusterData <- as.data.frame(colData(sce)[,c("Sample", clusterMethod)]) 
  
  # prepping the colors we want
  # for cell type
  num_pal <- length(unique(clusterData[,clusterMethod]))
  col_pal_ct <- grabColors(num_pal, start = 4)
  names(col_pal_ct) = unique(clusterData[,clusterMethod])
  # copying ct color pallete for Sample
  sample_pal <- length(unique(clusterData$Sample))
  col_pal_sample <- grabColors(sample_pal, start = 9)
  names(col_pal_sample) = unique(clusterData$Sample)
  
  # heatmap row annotationn
  row_ha <- rowAnnotation(
    df = clusterData[,c("Sample", clusterMethod)],
    col = list(cell_type = col_pal_ct)
  )
  
  heatmapped <- Heatmap(marker_z_score,
                        cluster_rows = TRUE,
                        cluster_columns = FALSE,
                        right_annotation = row_ha,
                        top_annotation = column_ha,
                        column_split = markTable$cellType,
                        column_title_rot = 30,
                        heatmap_legend_param = list(legend_gp = gpar(fontsize = 13)))
  
  print(heatmapped)
  
}

# printing 
pdf(here(plot_dir, "Markers_Heatmap_snAnno3_New_Pseudobulk.pdf"), width = 18, height = 22)
  pseudoHeater(sce, "snAnno3", markerList = progress.markers)
dev.off()

pdf(here(plot_dir, "Markers_Heatmap_snAnno3_New_Pseudobulk_SPLIT.pdf"), width = 18, height = 22)
  pseudoHeater(sce, "snAnno3", markerList = progress.markers_split)
dev.off()

pdf(here(plot_dir, "Markers_Heatmap_snAnno3_New_Pseudobulk_SPLIT_SINGLEgenes.pdf"), width = 18, height = 22)
  pseudoHeater(sce, "snAnno3", markerList = single.markers_split)
dev.off()

pdf(here(plot_dir, "Markers_Heatmap_snAnno3_New_Pseudobulk_SPLIT_SINGLEgenes_simplified.pdf"), width = 18, height = 22)
  pseudoHeater(sce, "snAnno3", markerList = single.markers_split_simplified)
dev.off()
# Done.