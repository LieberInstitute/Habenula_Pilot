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
library("RColorBrewer")

# creating plot directory
plot_dir <- here("plots", "misc_RPPR")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

load(here("processed-data", "04_snRNA-seq",  "sce_objects", "sce_pseudobulk_final_Annotations.Rdata"), verbose = TRUE)

# list of marker genes
official_markers = list(
  "Oligo" = c("MOBP"),
  "OPC" = c("PDGFRA"),
  "Micro" = c("CSF1R"),
  "Astro" = c("AQP4"),
  "Endo" = c("ITIH5"),
  "Thal" = c("LYPD6B", "ADARB2", "RORB"),
  # "LHb.1" = c("LINC02653"), #  , ATP8B1
  # "LHb.2" = c("AC073071.1"),
  # "LHb.3" = c ("ENTHD1"),
  # "LHb.4" = c("TLE2"),
  # "LHb.5" = c("LINC01619"),
  # "LHb.6" = c("TACR3"),
  # "LHb.7" = c("AC008619.1"),
  # "MHb.1" = c("EXOC1L"),
  # "MHb.2" = c("CHAT"),
  # "MHb.3" = c("BHLHE22"),
  'Hb' = c("POU4F1", "GPR151", "TAC3"), # BARHL1
  "MHb" = c("CHRNB4"),
  "LHb" = c("HTR2C"),
  'Neuron' = c('SYT1'),
  'Exc_Neuron' = c('SLC17A6'),
  'Inh_Neuron' = c('GAD1')
)

#### Check marker genes ####
marker_stats <- readxl::read_xlsx(here("plots", "99_paper_figs", "10c_snResolution_Top_Markers", "snResolution_top50MarkerGenes.xlsx"))
marker_stats[[1]] <- NULL
marker_stats |> count(cellType.target)

official_markers_tb <- tibble(cellType.short = names(official_markers),
                              cellType.target = names(official_markers),
                              gene = official_markers) |>
  unnest(gene) |>
  mutate(cellType.target = case_when(cellType.target == "Astro" ~'Astrocyte',
                                     cellType.target == "Micro" ~'Microglia',
                                     TRUE ~ cellType.target))

## Hb an Thal sub-type genes are datadriven (top mean ratio genes, glia is a mixed bag)
marker_details <- marker_stats |>
  right_join(official_markers_tb) |>
  arrange(rank_ratio) |>
  select(cellType.short, cellType.target, gene, rank_ratio, rank_marker) |>
  mutate(final_cell_type = cellType.target %in% marker_stats$cellType.target,
         anno = case_when(rank_ratio == 1 ~ "Data-Driven",
                          !final_cell_type | cellType.target == "Microglia" ~ "Literature",
                          TRUE ~ "Literature + Data-Supported"))

marker_details |> print(n = 22)

marker_stats |>
  group_by(cellType.target) |>
  arrange(rank_ratio) |>
  slice(1:5)

#### prep complex heatmap annotations  ####
# explicit color scheme
color_official_markers = c(
  "Oligo" = c("#4d5802"),
  "OPC"= c("#d3c871"),
  "Micro" = c("#222222"),
  "Astro" = c("#8d363c"),
  "Endo" = c("#ee6c14"),
  "Thal" = c("#EADDCA"),
  "LHb.1" = c("#0085af"),
  "LHb.2" = c("#0096FF"),
  "LHb.3" = c ("#89CFF0"),
  "LHb.4" = c("#6F8FAF"),
  "LHb.5" = c("#40E0D0"),
  "LHb.6" = c("#008080"),
  "LHb.7" = c("#7DF9FF"),
  "MHb.1" = c("#FF00FF"),
  "MHb.2" = c("#FAA0A0"),
  "MHb.3" = c("#fa246a"),
  "Hb" = c("#702963"),
  "MHb" = c("#F33A6A"),
  "LHb" = c("#0000FF"),
  'Neuron' = c("#5C5C5C"),
  'Exc_Neuron' = c('#8F8F8F'),
  "Inh_Neuron" = c("#C2C2C2")
)

marker_method_colors <- c(`Data-Driven` = "#FFDA85",
                          `Literature + Data-Supported` = "#F7A5A1",
                          `Literature` = "#95E4EE")

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
sce = sce_pb
clusterMethod = "final_Annotations"
markerList = official_markers

# Replacing genes with symbols for heatmap (remember, this is pseudobulked data)
rownames(sce) <- rowData(sce)$Symbol

# renaming rownnames of colData(sce) based on row annotations
rownames(colData(sce)) <- colData(sce)[, clusterMethod]

# Making data frame of genes we are interested in annd their general classification
# markTable <- as.data.frame(unlist(markerList)) |>
#   rownames_to_column("cellType") |>
#   rename(gene = `unlist(markerList)`) |>
#   # mutate(cellType = gsub("\\d+", "", cellType)) |>
#   filter(gene %in% rowData(sce_reorder)$Symbol)

markTable <- marker_details |>
    mutate(
        cellType.short = factor(cellType.short, levels = names(official_markers)),
        anno = factor(
            anno,
            levels = c("Data-Driven", "Literature + Data-Supported", "Literature")
        )
    ) |>
    arrange(anno, cellType.short)
## fix order in colsplit
official_markers_order <- c(
    "Hb",
    "MHb",
    "LHb",
    "Thal",
    "Neuron",
    "Exc_Neuron",
    "Inh_Neuron",
    "Oligo",
    "OPC",
    "Astro",
    "Endo",
    "Micro"
)
stopifnot(all(official_markers_order %in% names(official_markers)))
markTable$cellType.short = factor(
    markTable$cellType.short,
    levels = official_markers_order
)

markTable |> print(n=22)

row_namers <- c("LHb.1",
                "LHb.2",
                "LHb.3",
                "LHb.4",
                "LHb.5",
                "LHb.6",
                "LHb.7",
                "MHb.1",
                "MHb.2",
                "MHb.3",
                'Inhib.Thal',
                'Excit.Thal',
                "Oligo",
                "OPC",
                "Astrocyte",
                "Endo",
                "Microglia")

# row_namers <- markTable |> filter(final_cell_type) |> pull(cellType.target)
# marker_order <- markTable |> pull(gene)

# reordering sce object for plottability
sce_reorder <- sce[markTable$gene , row_namers]
sce_reorder$final_Annotations

# getting z scores
marker_z_score <- scale(t(logcounts(sce_reorder)))
# corner(marker_z_score)

identical(markTable$gene, colnames(marker_z_score))

# heatmap columns annotation
column_ha <- HeatmapAnnotation(
  Marker_Gene = markTable$cellType.short,
  # Marker_Method = markTable$anno,
  col = list(Marker_Gene = color_official_markers,
             Marker_Method = marker_method_colors)
)

# grabbing the annotations per cluster from the sce_reorder object
# clusterData <- as.data.frame(colData(sce_reorder)[,clusterMethod])
# names(clusterData) <- "cellType"

# prepping the colors we want
# for cell type
# num_pal <- length(unique(clusterData$cellType))
# col_pal_ct <- grabColors(num_pal, start = 4)
# names(col_pal_ct) = unique(clusterData$cellType)

# heatmap row annotationn
# identical(clusterData$cellType, rownames(marker_z_score))

row_ha <- rowAnnotation(
  Clusters = rownames(marker_z_score),
  col = list(Clusters = sn_colors)
)

heatmapped <- Heatmap(marker_z_score,
                      name = "Z Score",
                      col = circlize::colorRamp2(seq(-4, 4,  8/10),
                                                 rev(RColorBrewer::brewer.pal(11, "RdBu"))),
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      right_annotation = row_ha,
                      top_annotation = column_ha,
                      column_split = markTable$cellType.short,
                      column_title_rot = 30
                      # heatmap_legend_param = list(
                      #   title = c("Z_Score"),
                      #   border = "black"
                      # )
                      )

# printing
pdf(here(plot_dir, "Completed_Markers_Heatmap_final_Anno_FINAL.pdf"), width = 12, height = 8)
  heatmapped
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ──────────────────────────────────────────────────────────────────────────────
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
# ─ Packages ──────────────────────────────────────────────────────────────────────────────────
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
#
# R ── Package was removed from disk.
#
# ─────────────────────────────────────────────────────────────────────────────────────────────
