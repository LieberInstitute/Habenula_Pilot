## March 10, 2023 - Bukola Ajanaku
# Heatmapping updated pseudobulks with new marker lists
# qrsh -l mem_free=20G,h_vmem=20G
# Terminal 7

library("SingleCellExperiment")
library("here")
library("sessioninfo")
library("jaffelab")
library("ggplot2")
library("xlsx")
library("tidyverse")
library("ComplexHeatmap")

# grabbing new pseudobulked data 
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
            "sce_new_pseudobulk_with_snAnno_version2and3.Rdata"), verbose = TRUE)

# grabbing color schemes: grabColors(), max colors 20
source(here("code", "04_snRNA-seq", "sourcing", "color_Scheme_CT.R"))

# creating dir for plots
plot_dir <- here("plots", "05_explore_sce", "03_Plot_Investigate_snAnno", 
                 "03b_New_Pseudobulk_Heatmaps")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

# new markers custom
new_markers.custom <- list(
  'Neuron' = c('SYT1', 'SNAP25', "SYT4", "SYP"), 
  "Non-Neuronal Subtype" = c("TNF", "KIR4.1", "KCNJ10"),
  'excitatory_neuron' = c('SLC17A6', 'SLC17A7'), 
  'inhibitory_neuron' = c('GAD1', 'GAD2'), 
  "Habenula_neurons" = c("CCBP2","CD63", "HTR5B", "KCNNH8", "KCTD8", "LRRC55", "MAPK4", "NEUROD1", 
                         "PIXNC1", "SCUBE1", "SSTR4", "TACR1", "SSTR2", "IRX2"),
  "LHB_neuron_specific" = c("PCDH10", "GABRA1", "SYN2", "GAP43", "HTR2C", "ADCYAP1R1", 
                            "CHRM3", "VGF", "GPR151", "SST", "SLC17AR"),
  "MHb_neuron_specific" = c("TAC3", "SPON1", "SEMA3d", "CALB1",'HHTR5b', "CNR1", "GPR4", "CHRNB3", "B4",
                            "SLC17A7"),
  'mediodorsal thalamus'= c('EPHA4','PDYN', 'LYPD6B', 'LYPD6', 'S1PR1', 'GBX2', 
                            'RAMP3', 'COX6A2', 'SLITRK6', 'DGAT2', "ADARB2"),
  "Endo/Mural" = c("CLDN5", "CARMN", "ITIH5", "NOTCH3", "ATP10A", "MECOM", "EBF1", 
                   "AC092957.1", "ITGA1", "VWF"),
  'oligodendrocyte' = c('MOBP', 'MBP', "CX3CR1"),
  'oligodendrocyte_precursor' = c('PDGFRA', 'VCAN'), 
  'microglia' = c('C3', 'CSF1R'), 
  'astrocyte' = c('GFAP', 'AQP4')
)

### Creating function for heatmap with annotations
pseudoHeater <- function(sce, clusterMethod){
  # sce = pseudobulked sce object
  # clusterMethod = whichever nearest neighbors method you use (as character string)
  
  # Replacing genes with symbols for heatmap (remember, this is pseudobulked data)
  rownames(sce) <- rowData(sce)$Symbol
  
  # renaming rownnames of colData(sce) based on row annotations
  rownames(colData(sce)) <- paste(colData(sce)[, clusterMethod])
  
  # Making data frame of genes we are interested in annd their general classification
  markTable <- as.data.frame(unlist(new_markers.custom)) |> 
    rownames_to_column("cellType") |>
    rename(gene = `unlist(new_markers.custom)`) |>
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

# Running function and plotting 
## for snAnno2 which just has LHb.6 renames as Endo 
pdf(here(plot_dir, "Markers_Heatmap_snAnno2_Pseudobulked.pdf"), width = 18, height = 22)
  pseudoHeater(sce_new_pb_snAnno2, "snAnno2")
dev.off()

## for snAnno3 which has MHb.3 and MHb.2 combines
pdf(here(plot_dir, "Markers_Heatmap_snAnno3_New_Pseudobulk.pdf"), width = 18, height = 22)
  pseudoHeater(sce_new_pb_snAnno3, "snAnno3")
dev.off()



