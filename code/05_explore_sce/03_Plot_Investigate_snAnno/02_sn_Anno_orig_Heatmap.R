## March 8, 2023 - Bukola Ajanaku
# Taking the originally combined clusters from walktrap10 with the resolution 
# necessary for signle nucleus annotation in the lateral and medial habenula and 
# violin plotting against different sets of markers in order to clearly validate their 
# their identities.
# qrsh -l mem_free=20G,h_vmem=20G

library("SingleCellExperiment")
library("here")
library("sessioninfo")
library("jaffelab")
library("ggplot2")
library("xlsx")
library("tidyverse")
library("ComplexHeatmap")

####### LOAD/SOURCE #############################################################
# because this is for heatmps, this must be pseudobulked data (only need sce_psuedo_wT10)
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
          "sce_pseduobulked_wT.Rdata"), verbose = TRUE)

sce_pb <- sce_psuedo_wT10
sce_psuedo_wT20 <- NULL
sce_psuedo_wT50 <- NULL
  
# grabbing color schemes: grabColors(), max colors 20
source(here("code", "04_snRNA-seq", "sourcing", "color_Scheme_CT.R"))

# pulling excel sheet with annotations
annoWT10 <- read.xlsx(file = here("processed-data", "04_snRNA-seq", "09_Clustered_QC",
                                  "GeneAnnotations_Bukola.xlsx"), sheetName = c("WalkTrap10"))

# cleaning up for standardization
# 10
annoWT10_clean <- annoWT10 |> 
  filter(!is.na(Cluster)) |>
  mutate(snType_clean = gsub("[^a-zA-Z0-9]", "", snType), 
         Cluster = paste0("10wTrap_", Cluster))

# sanity check (looking for no repeat names in summarized list)
annoWT10_clean |> count(snType_clean)
    # 1     Astrocyte  3
    # 2     ExcitThal 11
    # 3            Hb  1
    # 4     InhibThal  5
    # 5          LHb1  1
    # 6          LHb2  1
    # 7          LHb3  1
    # 8          LHb4  1
    # 9          LHb5  1
    # 10         LHb6  1
    # 11         LHb7  1
    # 12         LHb8  1
    # 13         MHb1  1
    # 14         MHb2  1
    # 15         MHb3  1
    # 16         MHb4  1
    # 17    Microglia  1
    # 18          OPC  1
    # 19        Oligo  3

####### BODY ###################################################################
# using match to add annotated name columns for wt50
sce_pb$snAnno <- 
  annoWT10_clean$snType[match(sce_pb$wT_10_Erik, annoWT10_clean$Cluster)]

# list of marker genes 
markers.custom = list(
  'neuron' = c('SYT1', 'SNAP25'), 
  'excitatory_neuron' = c('SLC17A6', 'SLC17A7'), 
  'inhibitory_neuron' = c('GAD1', 'GAD2'), 
  'mediodorsal thalamus'= c('EPHA4','PDYN', 'LYPD6B', 'LYPD6', 'S1PR1', 'GBX2', 
                            'RAMP3', 'COX6A2', 'SLITRK6', 'DGAT2', "ADARB2"), 
  'Pf/PVT' = c("INHBA", "NPTXR"), 
  'Hb neuron specific'= c('POU4F1','GPR151', 'VAV2', "NR4A2", "NEUROD1"), 
  'MHB neuron specific' = c('TAC1','CHAT','CHRNB4', "TAC3"),
  'LHB neuron specific' = c('HTR2C','MMRN1', "ANO3", "ESRRG"),
  'oligodendrocyte' = c('MOBP', 'MBP'),
  'oligodendrocyte_precursor' = c('PDGFRA', 'VCAN'), 
  'microglia' = c('C3', 'CSF1R'), 
  'astrocyte' = c('GFAP', 'AQP4'),
  "Endo/Mural" = c("CLDN5", "CARMN", "ITIH5", "NOTCH3", "ATP10A", "MECOM", "EBF1", 
                   "AC092957.1", "ITGA1", "VWF"),
  "Choroid Plexus" = c("klotho", "CLIC6", "OATP14", "Ezrin")
)

####### PLOTTING ###############################################################
# creating dir for plots
plot_dir <- here("plots", "05_explore_sce", "03_Plot_Investigate_snAnno", 
                 "02_snAnno_orig_Heatmap")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

### Creating function for heatmap with annotations
pseudoHeater <- function(sce, clusterMethod){
  # sce = pseudobulked sce object
  # clusterMethod = whichever nearest neighbors method you use (as character string)
  
  # Replacing genes with symbols for heatmap (remember, this is pseudobulked data)
  rownames(sce) <- rowData(sce)$Symbol
  
  # renaming rownnames of colData(sce) based on row annotations
  rownames(colData(sce)) <- paste(colData(sce)[, "snAnno"], 
                                  ss(rownames(colData(sce)), "_", 3), sep = "_")
  
  # Making data frame of genes we are interested in annd their general classification
  markTable <- as.data.frame(unlist(markers.custom)) |> 
    rownames_to_column("cellType") |>
    rename(gene = `unlist(markers.custom)`) |>
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
    df = clusterData[,c("Sample", "snAnno")],
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

# Running function
pdf(here(plot_dir, "markers_heatmap_snAnno_original_pseudobulk_sce_obj.pdf"), width = 18, height = 22)
  pseudoHeater(sce_pb, "snAnno")
dev.off()


##### Using new pseudobulked sce object ########################################
## this is new because it pseudobulks our sce object based on my original snAnno 
# classificiation (before combining MHb 2&3). I will then heatmap to see the expression
# against my new custom markers list. 
# Terminal 7
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
          "sce_new_pseudobulk_with_snAnno.Rdata"))
new_sce_pb <- sce_new_pb_snAnno
sce_new_pb_snAnno <- NULL

# modified pseudoHeater for simplified pseudobulk
pseudoHeater <- function(sce, clusterMethod){
  # sce = pseudobulked sce object
  # clusterMethod = whichever nearest neighbors method you use (as character string)
  
  # Replacing genes with symbols for heatmap (remember, this is pseudobulked data)
  rownames(sce) <- rowData(sce)$Symbol
  
  # renaming rownnames of colData(sce) based on row annotations
  rownames(colData(sce)) <- paste(colData(sce)[, "snAnno"])
  
  # Making data frame of genes we are interested in annd their general classification
  markTable <- as.data.frame(unlist(markers.custom)) |> 
    rownames_to_column("cellType") |>
    rename(gene = `unlist(markers.custom)`) |>
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
    df = clusterData[,c("Sample", "snAnno")],
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

pdf(here(plot_dir, "markers_heatmap_snAnno_new_pseudobulk_sce_obj.pdf"), width = 18, height = 22)
  pseudoHeater(sce_pb, "snAnno")
dev.off()





# Done.

