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

## Reproducibility information
print("Reproducibility information:")
options(width = 120)
session_info()

# ─ Session info ─────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.2.3 Patched (2023-04-07 r84211)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2023-06-13
# pandoc   2.19.2 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/bin/pandoc
# 
# ─ Packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# Biobase              * 2.58.0    2022-11-01 [2] Bioconductor
# BiocGenerics         * 0.44.0    2022-11-01 [2] Bioconductor
# bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.2.1)
# circlize               0.4.15    2022-05-10 [2] CRAN (R 4.2.1)
# cli                    3.6.1     2023-03-23 [2] CRAN (R 4.2.3)
# clue                   0.3-64    2023-01-31 [2] CRAN (R 4.2.2)
# cluster                2.1.4     2022-08-22 [3] CRAN (R 4.2.3)
# codetools              0.2-19    2023-02-01 [3] CRAN (R 4.2.3)
# colorout             * 1.2-2     2023-02-09 [1] Github (jalvesaq/colorout@79931fd)
# colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.2.2)
# ComplexHeatmap       * 2.14.0    2022-11-01 [2] Bioconductor
# crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.2.1)
# DelayedArray           0.24.0    2022-11-01 [2] Bioconductor
# digest                 0.6.31    2022-12-11 [2] CRAN (R 4.2.2)
# doParallel             1.0.17    2022-02-07 [2] CRAN (R 4.2.1)
# dplyr                * 1.1.1     2023-03-22 [2] CRAN (R 4.2.3)
# fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.2.2)
# forcats              * 1.0.0     2023-01-29 [2] CRAN (R 4.2.2)
# foreach                1.5.2     2022-02-02 [2] CRAN (R 4.2.1)
# fs                     1.6.1     2023-02-06 [2] CRAN (R 4.2.2)
# gargle                 1.3.0     2023-01-30 [2] CRAN (R 4.2.2)
# generics               0.1.3     2022-07-05 [2] CRAN (R 4.2.1)
# GenomeInfoDb         * 1.34.9    2023-02-02 [2] Bioconductor
# GenomeInfoDbData       1.2.9     2022-09-29 [2] Bioconductor
# GenomicRanges        * 1.50.2    2022-12-16 [2] Bioconductor
# GetoptLong             1.0.5     2020-12-15 [2] CRAN (R 4.2.1)
# ggplot2              * 3.4.2     2023-04-03 [2] CRAN (R 4.2.3)
# GlobalOptions          0.1.2     2020-06-10 [2] CRAN (R 4.2.1)
# glue                   1.6.2     2022-02-24 [2] CRAN (R 4.2.1)
# googledrive            2.1.0     2023-03-22 [2] CRAN (R 4.2.3)
# gtable                 0.3.3     2023-03-21 [2] CRAN (R 4.2.3)
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.2.1)
# hms                    1.1.3     2023-03-21 [2] CRAN (R 4.2.3)
# IRanges              * 2.32.0    2022-11-01 [2] Bioconductor
# iterators              1.0.14    2022-02-05 [2] CRAN (R 4.2.1)
# jaffelab             * 0.99.32   2023-02-09 [1] Github (LieberInstitute/jaffelab@7b7afe3)
# lattice                0.20-45   2021-09-22 [3] CRAN (R 4.2.3)
# lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.2.1)
# limma                  3.54.2    2023-02-28 [2] Bioconductor
# lubridate            * 1.9.2     2023-02-10 [2] CRAN (R 4.2.2)
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.2.1)
# MASS                   7.3-58.2  2023-01-23 [3] CRAN (R 4.2.3)
# Matrix                 1.5-4     2023-04-04 [2] CRAN (R 4.2.3)
# MatrixGenerics       * 1.10.0    2022-11-01 [2] Bioconductor
# matrixStats          * 1.0.0     2023-06-02 [1] CRAN (R 4.2.3)
# munsell                0.5.0     2018-06-12 [2] CRAN (R 4.2.1)
# nlme                   3.1-162   2023-01-31 [2] CRAN (R 4.2.2)
# pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.2.3)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.2.1)
# png                    0.1-8     2022-11-29 [2] CRAN (R 4.2.2)
# purrr                * 1.0.1     2023-01-10 [2] CRAN (R 4.2.2)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.2.1)
# rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 4.2.2)
# RColorBrewer           1.1-3     2022-04-03 [2] CRAN (R 4.2.1)
# RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.2.3)
# readr                * 2.1.4     2023-02-10 [2] CRAN (R 4.2.2)
# rJava                  1.0-6     2021-12-10 [2] CRAN (R 4.2.1)
# rjson                  0.2.21    2022-01-09 [2] CRAN (R 4.2.1)
# rlang                  1.1.1     2023-04-28 [1] CRAN (R 4.2.3)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.2.1)
# S4Vectors            * 0.36.2    2023-02-26 [2] Bioconductor
# scales                 1.2.1     2022-08-20 [2] CRAN (R 4.2.1)
# segmented              1.6-4     2023-04-13 [1] CRAN (R 4.2.3)
# sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.2.1)
# shape                  1.4.6     2021-05-19 [2] CRAN (R 4.2.1)
# SingleCellExperiment * 1.20.1    2023-03-17 [2] Bioconductor
# stringi                1.7.12    2023-01-11 [2] CRAN (R 4.2.2)
# stringr              * 1.5.0     2022-12-02 [2] CRAN (R 4.2.2)
# SummarizedExperiment * 1.28.0    2022-11-01 [2] Bioconductor
# tibble               * 3.2.1     2023-03-20 [2] CRAN (R 4.2.3)
# tidyr                * 1.3.0     2023-01-24 [2] CRAN (R 4.2.2)
# tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.2.1)
# tidyverse            * 2.0.0     2023-02-22 [2] CRAN (R 4.2.2)
# timechange             0.2.0     2023-01-11 [2] CRAN (R 4.2.2)
# tzdb                   0.3.0     2022-03-28 [2] CRAN (R 4.2.1)
# utf8                   1.2.3     2023-01-31 [2] CRAN (R 4.2.2)
# vctrs                  0.6.2     2023-04-19 [1] CRAN (R 4.2.3)
# withr                  2.5.0     2022-03-03 [2] CRAN (R 4.2.1)
# xlsx                 * 0.6.5     2020-11-10 [2] CRAN (R 4.2.1)
# xlsxjars               0.6.1     2014-08-22 [2] CRAN (R 4.2.1)
# XVector                0.38.0    2022-11-01 [2] Bioconductor
# zlibbioc               1.44.0    2022-11-01 [2] Bioconductor
# 
# [1] /users/bsimbiat/R/4.2.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/library
