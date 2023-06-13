## 2/21/23 - Bukola Ajanaku
# Pseudobulking to on snAnno level to investigate new annotations 
# qrsh -l mem_free=20G,h_vmem=20G
# Terminal 7

library("SingleCellExperiment")
library("here")
library("spatialLIBD")
library("sessioninfo")

load(here("processed-data", "04_snRNA-seq", "sce_objects", 
          "sce_post_09_clustered_qc.Rdata"))

# updating LHb.6 as Endo 
# creating snAnno2 which will contain the combined MHb groups
sce$snAnno2 <- sce$snAnno

## LHb.6 is actually Endothelial. Total LHb is now 7 from 8.
sce$snAnno2[sce$snAnno2 == "LHb.6"] <- "Endo"

## Combining MHb.3 with MHb.2
sce$snAnno3 <- sce$snAnno2

# combining MHb.3 with MHb.2
sce$snAnno3[sce$snAnno3 == "MHb.3"] <- "MHb.2"


# Pseudobulking to create compressed sce object
message("Start - snAnno2")
Sys.time()

  sce_new_pb_snAnno2 <- registration_pseudobulk(sce, "snAnno2", "Sample")

message("End - Pseudobulk Complete")
Sys.time()

message("Start - snAnno3")
Sys.time()

  sce_new_pb_snAnno3 <- registration_pseudobulk(sce, "snAnno3", "Sample")

message("End - Pseudobulk Complete")
Sys.time()

save(sce_new_pb_snAnno2, sce_new_pb_snAnno3, file = here("processed-data", "04_snRNA-seq", "sce_objects", 
                                            "sce_new_pseudobulk_with_snAnno_version2and3.Rdata"))

# sgejobs::job_single('07c_Pseudobulking', create_shell = TRUE, queue= 'bluejay', memory = '25G', command = "Rscript 07c_Pseudobulking.R")

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
# cli                      3.6.1     2023-03-23 [2] CRAN (R 4.2.3)
# codetools                0.2-19    2023-02-01 [3] CRAN (R 4.2.3)
# colorout               * 1.2-2     2023-02-09 [1] Github (jalvesaq/colorout@79931fd)
# colorspace               2.1-0     2023-01-23 [2] CRAN (R 4.2.2)
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
# dplyr                    1.1.1     2023-03-22 [2] CRAN (R 4.2.3)
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
# foreach                  1.5.2     2022-02-02 [2] CRAN (R 4.2.1)
# generics                 0.1.3     2022-07-05 [2] CRAN (R 4.2.1)
# GenomeInfoDb           * 1.34.9    2023-02-02 [2] Bioconductor
# GenomeInfoDbData         1.2.9     2022-09-29 [2] Bioconductor
# GenomicAlignments        1.34.1    2023-03-09 [2] Bioconductor
# GenomicRanges          * 1.50.2    2022-12-16 [2] Bioconductor
# ggbeeswarm               0.7.1     2022-12-16 [2] CRAN (R 4.2.2)
# ggplot2                  3.4.2     2023-04-03 [2] CRAN (R 4.2.3)
# ggrepel                  0.9.3     2023-02-03 [2] CRAN (R 4.2.2)
# glue                     1.6.2     2022-02-24 [2] CRAN (R 4.2.1)
# golem                    0.4.0     2023-03-12 [2] CRAN (R 4.2.3)
# gridExtra                2.3       2017-09-09 [2] CRAN (R 4.2.1)
# gtable                   0.3.3     2023-03-21 [2] CRAN (R 4.2.3)
# HDF5Array                1.26.0    2022-11-01 [2] Bioconductor
# here                   * 1.0.1     2020-12-13 [2] CRAN (R 4.2.1)
# htmltools                0.5.5     2023-03-23 [2] CRAN (R 4.2.3)
# htmlwidgets              1.6.2     2023-03-17 [2] CRAN (R 4.2.3)
# httpuv                   1.6.9     2023-02-14 [2] CRAN (R 4.2.2)
# httr                     1.4.5     2023-02-24 [2] CRAN (R 4.2.2)
# interactiveDisplayBase   1.36.0    2022-11-01 [2] Bioconductor
# IRanges                * 2.32.0    2022-11-01 [2] Bioconductor
# irlba                    2.3.5.1   2022-10-03 [2] CRAN (R 4.2.1)
# iterators                1.0.14    2022-02-05 [2] CRAN (R 4.2.1)
# jquerylib                0.1.4     2021-04-26 [2] CRAN (R 4.2.1)
# jsonlite                 1.8.5     2023-06-05 [1] CRAN (R 4.2.3)
# KEGGREST                 1.38.0    2022-11-01 [2] Bioconductor
# later                    1.3.0     2021-08-18 [2] CRAN (R 4.2.1)
# lattice                  0.20-45   2021-09-22 [3] CRAN (R 4.2.3)
# lazyeval                 0.2.2     2019-03-15 [2] CRAN (R 4.2.1)
# lifecycle                1.0.3     2022-10-07 [2] CRAN (R 4.2.1)
# limma                    3.54.2    2023-02-28 [2] Bioconductor
# locfit                   1.5-9.7   2023-01-02 [2] CRAN (R 4.2.2)
# magick                   2.7.4     2023-03-09 [2] CRAN (R 4.2.3)
# magrittr                 2.0.3     2022-03-30 [2] CRAN (R 4.2.1)
# maps                     3.4.1     2022-10-30 [2] CRAN (R 4.2.2)
# Matrix                   1.5-4     2023-04-04 [2] CRAN (R 4.2.3)
# MatrixGenerics         * 1.10.0    2022-11-01 [2] Bioconductor
# matrixStats            * 1.0.0     2023-06-02 [1] CRAN (R 4.2.3)
# memoise                  2.0.1     2021-11-26 [2] CRAN (R 4.2.1)
# mime                     0.12      2021-09-28 [2] CRAN (R 4.2.1)
# munsell                  0.5.0     2018-06-12 [2] CRAN (R 4.2.1)
# paletteer                1.5.0     2022-10-19 [2] CRAN (R 4.2.1)
# pillar                   1.9.0     2023-03-22 [2] CRAN (R 4.2.3)
# pkgconfig                2.0.3     2019-09-22 [2] CRAN (R 4.2.1)
# plotly                   4.10.1    2022-11-07 [2] CRAN (R 4.2.2)
# png                      0.1-8     2022-11-29 [2] CRAN (R 4.2.2)
# promises                 1.2.0.1   2021-02-11 [2] CRAN (R 4.2.1)
# purrr                    1.0.1     2023-01-10 [2] CRAN (R 4.2.2)
# R.methodsS3              1.8.2     2022-06-13 [2] CRAN (R 4.2.1)
# R.oo                     1.25.0    2022-06-12 [2] CRAN (R 4.2.1)
# R.utils                  2.12.2    2022-11-11 [2] CRAN (R 4.2.2)
# R6                       2.5.1     2021-08-19 [2] CRAN (R 4.2.1)
# rappdirs                 0.3.3     2021-01-31 [2] CRAN (R 4.2.1)
# RColorBrewer             1.1-3     2022-04-03 [2] CRAN (R 4.2.1)
# Rcpp                     1.0.10    2023-01-22 [2] CRAN (R 4.2.2)
# RCurl                    1.98-1.12 2023-03-27 [2] CRAN (R 4.2.3)
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
# sessioninfo            * 1.2.2     2021-12-06 [2] CRAN (R 4.2.1)
# shiny                    1.7.4     2022-12-15 [2] CRAN (R 4.2.2)
# shinyWidgets             0.7.6     2023-01-08 [2] CRAN (R 4.2.2)
# SingleCellExperiment   * 1.20.1    2023-03-17 [2] Bioconductor
# spam                     2.9-1     2022-08-07 [2] CRAN (R 4.2.1)
# sparseMatrixStats        1.10.0    2022-11-01 [2] Bioconductor
# SpatialExperiment      * 1.8.1     2023-03-05 [2] Bioconductor
# spatialLIBD            * 1.10.1    2022-12-01 [2] Bioconductor
# statmod                  1.5.0     2023-01-06 [2] CRAN (R 4.2.2)
# SummarizedExperiment   * 1.28.0    2022-11-01 [2] Bioconductor
# tibble                   3.2.1     2023-03-20 [2] CRAN (R 4.2.3)
# tidyr                    1.3.0     2023-01-24 [2] CRAN (R 4.2.2)
# tidyselect               1.2.0     2022-10-10 [2] CRAN (R 4.2.1)
# utf8                     1.2.3     2023-01-31 [2] CRAN (R 4.2.2)
# vctrs                    0.6.2     2023-04-19 [1] CRAN (R 4.2.3)
# vipor                    0.4.5     2017-03-22 [2] CRAN (R 4.2.1)
# viridis                  0.6.2     2021-10-13 [2] CRAN (R 4.2.1)
# viridisLite              0.4.2     2023-05-02 [1] CRAN (R 4.2.3)
# XML                      3.99-0.14 2023-03-19 [2]

