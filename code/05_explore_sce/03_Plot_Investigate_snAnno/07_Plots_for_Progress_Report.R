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

# getting sce object
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
          "sce_post_09_clustered_qc.Rdata"))

# creating plot_dir
plot_dir <- here("plots", "05_explore_sce", "03_Plot_Investigate_snAnno", 
                 "07_Plots_for_Progress_Report")
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
color_row_namers <- c( "Oligo" = c("#5C4033"), # dark brown
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

  
# finalizing sce for ultimate save
sce_final <- sce
sce_final$final_Annotations <- sce$snAnno3

# printing 
pdf(here(plot_dir, "Completed_Markers_Heatmap_snAnno3_Simple_Pseudobulk.pdf"), width = 12, height = 8)
  print(heatmapped)
dev.off()

# saving official sce object
save(sce_final, file = here("processed-data", "04_snRNA-seq", "sce_objects", 
          "sce_final.Rdata"))


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
# ! package                * version   date (UTC) lib source
# AnnotationDbi            1.60.2    2023-03-10 [2] Bioconductor
# AnnotationHub            3.6.0     2022-11-01 [2] Bioconductor
# attempt                  0.3.1     2020-05-03 [2] CRAN (R 4.2.1)
# beachmat                 2.14.2    2023-04-07 [2] Bioconductor
# beeswarm                 0.4.0     2021-06-01 [2] CRAN (R 4.2.1)
# benchmarkme              1.0.8     2022-06-12 [2] CRAN (R 4.2.1)
# benchmarkmeData          1.0.4     2020-04-23 [2] CRAN (R 4.2.1)
# Biobase                * 2.58.0    2022-11-01 [2] Bioconductor
# BiocFileCache            2.6.1     2023-02-17 [2] Bioconductor
# BiocGenerics           * 0.44.0    2022-11-01 [2] Bioconductor
# BiocIO                   1.8.0     2022-11-01 [2] Bioconductor
# BiocManager              1.30.20   2023-02-24 [2] CRAN (R 4.2.2)
# BiocNeighbors            1.16.0    2022-11-01 [2] Bioconductor
# BiocParallel             1.32.6    2023-03-17 [2] Bioconductor
# BiocSingular             1.14.0    2022-11-01 [2] Bioconductor
# BiocVersion              3.16.0    2022-04-26 [2] Bioconductor
# Biostrings               2.66.0    2022-11-01 [2] Bioconductor
# bit                      4.0.5     2022-11-15 [2] CRAN (R 4.2.2)
# bit64                    4.0.5     2020-08-30 [2] CRAN (R 4.2.1)
# bitops                   1.0-7     2021-04-24 [2] CRAN (R 4.2.1)
# blob                     1.2.4     2023-03-17 [2] CRAN (R 4.2.3)
# bslib                    0.4.2     2022-12-16 [2] CRAN (R 4.2.2)
# cachem                   1.0.8     2023-05-01 [1] CRAN (R 4.2.3)
# circlize                 0.4.15    2022-05-10 [2] CRAN (R 4.2.1)
# cli                      3.6.1     2023-03-23 [2] CRAN (R 4.2.3)
# clue                     0.3-64    2023-01-31 [2] CRAN (R 4.2.2)
# cluster                  2.1.4     2022-08-22 [3] CRAN (R 4.2.3)
# codetools                0.2-19    2023-02-01 [3] CRAN (R 4.2.3)
# colorout               * 1.2-2     2023-02-09 [1] Github (jalvesaq/colorout@79931fd)
# colorspace               2.1-0     2023-01-23 [2] CRAN (R 4.2.2)
# ComplexHeatmap         * 2.14.0    2022-11-01 [2] Bioconductor
# config                   0.3.1     2020-12-17 [2] CRAN (R 4.2.1)
# cowplot                  1.1.1     2020-12-30 [2] CRAN (R 4.2.1)
# crayon                   1.5.2     2022-09-29 [2] CRAN (R 4.2.1)
# curl                     5.0.0     2023-01-12 [2] CRAN (R 4.2.2)
# data.table               1.14.8    2023-02-17 [2] CRAN (R 4.2.2)
# DBI                      1.1.3     2022-06-18 [2] CRAN (R 4.2.1)
# dbplyr                   2.3.2     2023-03-21 [2] CRAN (R 4.2.3)
# DelayedArray             0.24.0    2022-11-01 [2] Bioconductor
# DelayedMatrixStats       1.20.0    2022-11-01 [2] Bioconductor
# digest                   0.6.31    2022-12-11 [2] CRAN (R 4.2.2)
# doParallel               1.0.17    2022-02-07 [2] CRAN (R 4.2.1)
# dotCall64                1.0-2     2022-10-03 [2] CRAN (R 4.2.1)
# dplyr                  * 1.1.1     2023-03-22 [2] CRAN (R 4.2.3)
# dqrng                    0.3.0     2021-05-01 [2] CRAN (R 4.2.1)
# DropletUtils             1.18.1    2022-11-22 [2] Bioconductor
# DT                       0.27      2023-01-17 [2] CRAN (R 4.2.2)
# edgeR                    3.40.2    2023-01-19 [2] Bioconductor
# ellipsis                 0.3.2     2021-04-29 [2] CRAN (R 4.2.1)
# ExperimentHub            2.6.0     2022-11-01 [2] Bioconductor
# fansi                    1.0.4     2023-01-22 [2] CRAN (R 4.2.2)
# fastmap                  1.1.1     2023-02-24 [2] CRAN (R 4.2.2)
# fields                   14.1      2022-08-12 [2] CRAN (R 4.2.1)
# filelock                 1.0.2     2018-10-05 [2] CRAN (R 4.2.1)
# forcats                * 1.0.0     2023-01-29 [2] CRAN (R 4.2.2)
# foreach                  1.5.2     2022-02-02 [2] CRAN (R 4.2.1)
# fs                       1.6.1     2023-02-06 [2] CRAN (R 4.2.2)
# gargle                   1.3.0     2023-01-30 [2] CRAN (R 4.2.2)
# generics                 0.1.3     2022-07-05 [2] CRAN (R 4.2.1)
# GenomeInfoDb           * 1.34.9    2023-02-02 [2] Bioconductor
# GenomeInfoDbData         1.2.9     2022-09-29 [2] Bioconductor
# GenomicAlignments        1.34.1    2023-03-09 [2] Bioconductor
# GenomicRanges          * 1.50.2    2022-12-16 [2] Bioconductor
# GetoptLong               1.0.5     2020-12-15 [2] CRAN (R 4.2.1)
# ggbeeswarm               0.7.1     2022-12-16 [2] CRAN (R 4.2.2)
# ggplot2                * 3.4.2     2023-04-03 [2] CRAN (R 4.2.3)
# ggrepel                  0.9.3     2023-02-03 [2] CRAN (R 4.2.2)
# GlobalOptions            0.1.2     2020-06-10 [2] CRAN (R 4.2.1)
# glue                     1.6.2     2022-02-24 [2] CRAN (R 4.2.1)
# golem                    0.4.0     2023-03-12 [2] CRAN (R 4.2.3)
# googledrive              2.1.0     2023-03-22 [2] CRAN (R 4.2.3)
# gridExtra                2.3       2017-09-09 [2] CRAN (R 4.2.1)
# gtable                   0.3.3     2023-03-21 [2] CRAN (R 4.2.3)
# HDF5Array                1.26.0    2022-11-01 [2] Bioconductor
# here                   * 1.0.1     2020-12-13 [2] CRAN (R 4.2.1)
# hms                      1.1.3     2023-03-21 [2] CRAN (R 4.2.3)
# htmltools                0.5.5     2023-03-23 [2] CRAN (R 4.2.3)
# htmlwidgets              1.6.2     2023-03-17 [2] CRAN (R 4.2.3)
# httpuv                   1.6.9     2023-02-14 [2] CRAN (R 4.2.2)
# httr                     1.4.5     2023-02-24 [2] CRAN (R 4.2.2)
# interactiveDisplayBase   1.36.0    2022-11-01 [2] Bioconductor
# IRanges                * 2.32.0    2022-11-01 [2] Bioconductor
# irlba                    2.3.5.1   2022-10-03 [2] CRAN (R 4.2.1)
# iterators                1.0.14    2022-02-05 [2] CRAN (R 4.2.1)
# jaffelab               * 0.99.32   2023-02-09 [1] Github (LieberInstitute/jaffelab@7b7afe3)
# jquerylib                0.1.4     2021-04-26 [2] CRAN (R 4.2.1)
# jsonlite                 1.8.5     2023-06-05 [1] CRAN (R 4.2.3)
# KEGGREST                 1.38.0    2022-11-01 [2] Bioconductor
# later                    1.3.0     2021-08-18 [2] CRAN (R 4.2.1)
# lattice                  0.20-45   2021-09-22 [3] CRAN (R 4.2.3)
# lazyeval                 0.2.2     2019-03-15 [2] CRAN (R 4.2.1)
# lifecycle                1.0.3     2022-10-07 [2] CRAN (R 4.2.1)
# limma                    3.54.2    2023-02-28 [2] Bioconductor
# locfit                   1.5-9.7   2023-01-02 [2] CRAN (R 4.2.2)
# lubridate              * 1.9.2     2023-02-10 [2] CRAN (R 4.2.2)
# magick                   2.7.4     2023-03-09 [2] CRAN (R 4.2.3)
# magrittr                 2.0.3     2022-03-30 [2] CRAN (R 4.2.1)
# maps                     3.4.1     2022-10-30 [2] CRAN (R 4.2.2)
# MASS                     7.3-58.2  2023-01-23 [3] CRAN (R 4.2.3)
# Matrix                   1.5-4     2023-04-04 [2] CRAN (R 4.2.3)
# MatrixGenerics         * 1.10.0    2022-11-01 [2] Bioconductor
# matrixStats            * 1.0.0     2023-06-02 [1] CRAN (R 4.2.3)
# memoise                  2.0.1     2021-11-26 [2] CRAN (R 4.2.1)
# mime                     0.12      2021-09-28 [2] CRAN (R 4.2.1)
# munsell                  0.5.0     2018-06-12 [2] CRAN (R 4.2.1)
# nlme                     3.1-162   2023-01-31 [2] CRAN (R 4.2.2)
# paletteer                1.5.0     2022-10-19 [2] CRAN (R 4.2.1)
# pillar                   1.9.0     2023-03-22 [2] CRAN (R 4.2.3)
# pkgconfig                2.0.3     2019-09-22 [2] CRAN (R 4.2.1)
# plotly                   4.10.1    2022-11-07 [2] CRAN (R 4.2.2)
# png                      0.1-8     2022-11-29 [2] CRAN (R 4.2.2)
# promises                 1.2.0.1   2021-02-11 [2] CRAN (R 4.2.1)
# purrr                  * 1.0.1     2023-01-10 [2] CRAN (R 4.2.2)
# R.methodsS3              1.8.2     2022-06-13 [2] CRAN (R 4.2.1)
# R.oo                     1.25.0    2022-06-12 [2] CRAN (R 4.2.1)
# R.utils                  2.12.2    2022-11-11 [2] CRAN (R 4.2.2)
# R6                       2.5.1     2021-08-19 [2] CRAN (R 4.2.1)
# rafalib                * 1.0.0     2015-08-09 [1] CRAN (R 4.2.2)
# rappdirs                 0.3.3     2021-01-31 [2] CRAN (R 4.2.1)
# RColorBrewer             1.1-3     2022-04-03 [2] CRAN (R 4.2.1)
# Rcpp                     1.0.10    2023-01-22 [2] CRAN (R 4.2.2)
# RCurl                    1.98-1.12 2023-03-27 [2] CRAN (R 4.2.3)
# readr                  * 2.1.4     2023-02-10 [2] CRAN (R 4.2.2)
# rematch2                 2.1.2     2020-05-01 [2] CRAN (R 4.2.1)
# restfulr                 0.0.15    2022-06-16 [2] CRAN (R 4.2.1)
# rhdf5                    2.42.1    2023-04-07 [2] Bioconductor
# rhdf5filters             1.10.1    2023-03-24 [2] Bioconductor
# Rhdf5lib                 1.20.0    2022-11-01 [2] Bioconductor
# rjson                    0.2.21    2022-01-09 [2] CRAN (R 4.2.1)
# rlang                    1.1.1     2023-04-28 [1] CRAN (R 4.2.3)
# rprojroot                2.0.3     2022-04-02 [2] CRAN (R 4.2.1)
# Rsamtools                2.14.0    2022-11-01 [2] Bioconductor
# RSQLite                  2.3.1     2023-04-03 [2] CRAN (R 4.2.3)
# rsvd                     1.0.5     2021-04-16 [2] CRAN (R 4.2.1)
# rtracklayer              1.58.0    2022-11-01 [2] Bioconductor
# S4Vectors              * 0.36.2    2023-02-26 [2] Bioconductor
# R sass                     0.4.5     <NA>       [2] <NA>
#   ScaledMatrix             1.6.0     2022-11-01 [2] Bioconductor
# scales                   1.2.1     2022-08-20 [2] CRAN (R 4.2.1)
# scater                   1.26.1    2022-11-13 [2] Bioconductor
# scuttle                  1.8.4     2023-01-19 [2] Bioconductor
# segmented                1.6-4     2023-04-13 [1] CRAN (R 4.2.3)
# sessioninfo            * 1.2.2     2021-12-06 [2] CRAN (R 4.2.1)
# shape                    1.4.6     2021-05-19 [2] CRAN (R 4.2.1)
# shiny                    1.7.4     2022-12-15 [2] CRAN (R 4.2.2)
# shinyWidgets             0.7.6     2023-01-08 [2] CRAN (R 4.2.2)
# SingleCellExperiment   * 1.20.1    2023-03-17 [2] Bioconductor
# spam                     2.9-1     2022-08-07 [2] CRAN (R 4.2.1)
# sparseMatrixStats        1.10.0    2022-11-01 [2] Bioconductor
# SpatialExperiment      * 1.8.1     2023-03-05 [2] Bioconductor
# spatialLIBD            * 1.10.1    2022-12-01 [2] Bioconductor
# statmod                  1.5.0     2023-01-06 [2] CRAN (R 4.2.2)
# stringi                  1.7.12    2023-01-11 [2] CRAN (R 4.2.2)
# stringr                * 1.5.0     2022-12-02 [2] CRAN (R 4.2.2)
# SummarizedExperiment   * 1.28.0    2022-11-01 [2] Bioconductor
# tibble                 * 3.2.1     2023-03-20 [2] CRAN (R 4.2.3)
# tidyr                  * 1.3.0     2023-01-24 [2] CRAN (R 4.2.2)
# tidyselect               1.2.0     2022-10-10 [2] CRAN (R 4.2.1)
# tidyverse              * 2.0.0     2023-02-22 [2] CRAN (R 4.2.2)
# timechange               0.2.0     2023-01-11 [2] CRAN (R 4.2.2)
# tzdb                     0.3.0     2022-03-28 [2] CRAN (R 4.2.1)
# utf8                     1.2.3     2023-01-31 [2] CRAN (R 4.2.2)
# vctrs                    0.6.2     2023-04-19 [1] CRAN (R 4.2.3)
# vipor                    0.4.5     2017-03-22 [2] CRAN (R 4.2.1)
# viridis                  0.6.2     2021-10-13 [2] CRAN (R 4.2.1)
# viridisLite              0.4.2     2023-05-02 [1] CRAN (R 4.2.3)
# withr                    2.5.0     2022-03-03 [2] CRAN (R 4.2.1)
# XML                      3.99-0.14 2023-03-19 [2] CRAN (R 4.2.3)
# xtable                   1.8-4     2019-04-21 [2] CRAN (R 4.2.1)
# XVector                  0.38.0    2022-11-01 [2] Bioconductor
# yaml                     2.3.7     2023-01-23 [2] CRAN (R 4.2.2)
# zlibbioc                 1.44.0    2022-11-01 [2] Bioconductor
# 
# [1] /users/bsimbiat/R/4.2.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/library



