## June 11, 2023 - Bukola Ajanaku
# Creating official Wallace et al. sce object with cluster info and annotations 
# because they were not in the original sce object and were really hard to find.
# qrsh -l mem_free=50G,h_vmem=50G

library("SingleCellExperiment")
library("here")
library("sessioninfo")
library("jaffelab")
library("Seurat")
library("scater")
library("EnsDb.Hsapiens.v86")
#BiocManager::install("org.Rn.eg.db")
library("org.Hs.eg.db")
# library("RCurl")
library("AnnotationHub")
library("xlsx")

# loading Wallace et al. clean data, updating to version 3 of seurat object so 
# I can just make into a sce object 
wallData <- as.SingleCellExperiment(
  UpdateSeuratObject(
    readRDS(file = here("processed-data", 
                        "09_trans_special_analysis",
                        "Wallace_mouse_data.rds"))
  ))

wallSeu <- readRDS(file = here("processed-data", 
                               "09_trans_special_analysis",
                               "Wallace_mouse_data.rds"))

wall2s <- UpdateSeuratObject(wallSeu)
idents(wallSeu)

############ ADDING SAMPLE INFORMATION #########################################
# realized that the sample data was not in the colData of the Wallace sce object
# last minute
colData(wallData)$Sample <- NA

# grabbing nuclei identifiers
colData(wallData)$Row <- rownames(colData(wallData))

# adding Sample info
colData(wallData)[startsWith(colData(wallData)$Row, "hab_160822" ), ]$Sample <- "Mouse1"
colData(wallData)[startsWith(colData(wallData)$Row, "hab_161102" ), ]$Sample <- "Mouse2"
colData(wallData)[startsWith(colData(wallData)$Row, "hab_161103" ), ]$Sample <- "Mouse3"
colData(wallData)[startsWith(colData(wallData)$Row, "hab_161105" ), ]$Sample <- "Mouse4"

# adding hemispheric information
colData(wallData)$Hemi <- NA
colData(wallData)$Hemi <- ss(x = colData(wallData)$Row, "_", slot = 4)

############ ADDING CLUSTER INFORMATION #########################################
# adding Wallace cell types (found in supplements of eLife paper.)
eLifeCT <- read.csv(file = here("processed-data", "09_trans_special_analysis",
                                "elife_ct.csv"))

# grabbing cluster info
eLife_clusts <- read.xlsx(file = here("processed-data", "09_trans_special_analysis",
                                "elife_clusters.xlsx"), 1, header=TRUE)

names(eLife_clusts)[[1]] <- "gene"

# adding gene names to cell type data
# grabbing Ensembl GRCm38 release 87 information for Wallace data set
ah <- AnnotationHub()
query(ah, "EnsDb.Mmusculus.v87")
org.Mm.eg.db <- ah[["AH53222"]]

# grabbing relevant data from ensdb
mouse_genes <- DataFrame(genes(org.Mm.eg.db))

# making sure wallData has a column for it's symbols
eLifeCT$gene <- rownames(wallData)

# grabbing symbol IDs
eLifeCT$gene_id <- mouse_genes$gene_id[match(eLifeCT$gene, mouse_genes$gene_name)]

table(is.na(eLifeCT$gene_id))
# FALSE  TRUE 
# 22944  2345

## Reproducibility information
print("Reproducibility information:")
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
# package                * version   date (UTC) lib source
# abind                    1.4-5     2016-07-21 [2] CRAN (R 4.2.1)
# AnnotationDbi          * 1.60.2    2023-03-10 [2] Bioconductor
# AnnotationFilter       * 1.22.0    2022-11-01 [2] Bioconductor
# AnnotationHub          * 3.6.0     2022-11-01 [2] Bioconductor
# beachmat                 2.14.2    2023-04-07 [2] Bioconductor
# beeswarm                 0.4.0     2021-06-01 [2] CRAN (R 4.2.1)
# Biobase                * 2.58.0    2022-11-01 [2] Bioconductor
# BiocFileCache          * 2.6.1     2023-02-17 [2] Bioconductor
# BiocGenerics           * 0.44.0    2022-11-01 [2] Bioconductor
# BiocIO                   1.8.0     2022-11-01 [2] Bioconductor
# BiocManager              1.30.20   2023-02-24 [2] CRAN (R 4.2.2)
# BiocNeighbors            1.16.0    2022-11-01 [2] Bioconductor
# BiocParallel             1.32.6    2023-03-17 [2] Bioconductor
# BiocSingular             1.14.0    2022-11-01 [2] Bioconductor
# BiocVersion              3.16.0    2022-04-26 [2] Bioconductor
# biomaRt                  2.54.1    2023-03-20 [2] Bioconductor
# Biostrings               2.66.0    2022-11-01 [2] Bioconductor
# bit                      4.0.5     2022-11-15 [2] CRAN (R 4.2.2)
# bit64                    4.0.5     2020-08-30 [2] CRAN (R 4.2.1)
# bitops                   1.0-7     2021-04-24 [2] CRAN (R 4.2.1)
# blob                     1.2.4     2023-03-17 [2] CRAN (R 4.2.3)
# cachem                   1.0.8     2023-05-01 [1] CRAN (R 4.2.3)
# cli                      3.6.1     2023-03-23 [2] CRAN (R 4.2.3)
# cluster                  2.1.4     2022-08-22 [3] CRAN (R 4.2.3)
# codetools                0.2-19    2023-02-01 [3] CRAN (R 4.2.3)
# colorout               * 1.2-2     2023-02-09 [1] Github (jalvesaq/colorout@79931fd)
# colorspace               2.1-0     2023-01-23 [2] CRAN (R 4.2.2)
# cowplot                  1.1.1     2020-12-30 [2] CRAN (R 4.2.1)
# crayon                   1.5.2     2022-09-29 [2] CRAN (R 4.2.1)
# curl                     5.0.0     2023-01-12 [2] CRAN (R 4.2.2)
# data.table               1.14.8    2023-02-17 [2] CRAN (R 4.2.2)
# DBI                      1.1.3     2022-06-18 [2] CRAN (R 4.2.1)
# dbplyr                 * 2.3.2     2023-03-21 [2] CRAN (R 4.2.3)
# DelayedArray             0.24.0    2022-11-01 [2] Bioconductor
# DelayedMatrixStats       1.20.0    2022-11-01 [2] Bioconductor
# deldir                   1.0-6     2021-10-23 [2] CRAN (R 4.2.1)
# digest                   0.6.31    2022-12-11 [2] CRAN (R 4.2.2)
# dplyr                    1.1.1     2023-03-22 [2] CRAN (R 4.2.3)
# ellipsis                 0.3.2     2021-04-29 [2] CRAN (R 4.2.1)
# EnsDb.Hsapiens.v86     * 2.99.0    2023-02-09 [1] Bioconductor
# ensembldb              * 2.22.0    2022-11-01 [2] Bioconductor
# fansi                    1.0.4     2023-01-22 [2] CRAN (R 4.2.2)
# fastmap                  1.1.1     2023-02-24 [2] CRAN (R 4.2.2)
# filelock                 1.0.2     2018-10-05 [2] CRAN (R 4.2.1)
# fitdistrplus             1.1-11    2023-04-25 [1] CRAN (R 4.2.3)
# fs                       1.6.1     2023-02-06 [2] CRAN (R 4.2.2)
# future                   1.32.0    2023-03-07 [2] CRAN (R 4.2.3)
# future.apply             1.11.0    2023-05-21 [1] CRAN (R 4.2.3)
# gargle                   1.3.0     2023-01-30 [2] CRAN (R 4.2.2)
# generics                 0.1.3     2022-07-05 [2] CRAN (R 4.2.1)
# GenomeInfoDb           * 1.34.9    2023-02-02 [2] Bioconductor
# GenomeInfoDbData         1.2.9     2022-09-29 [2] Bioconductor
# GenomicAlignments        1.34.1    2023-03-09 [2] Bioconductor
# GenomicFeatures        * 1.50.4    2023-01-24 [2] Bioconductor
# GenomicRanges          * 1.50.2    2022-12-16 [2] Bioconductor
# ggbeeswarm               0.7.1     2022-12-16 [2] CRAN (R 4.2.2)
# ggplot2                * 3.4.2     2023-04-03 [2] CRAN (R 4.2.3)
# ggrepel                  0.9.3     2023-02-03 [2] CRAN (R 4.2.2)
# ggridges                 0.5.4     2022-09-26 [2] CRAN (R 4.2.1)
# globals                  0.16.2    2022-11-21 [2] CRAN (R 4.2.2)
# glue                     1.6.2     2022-02-24 [2] CRAN (R 4.2.1)
# goftest                  1.2-3     2021-10-07 [1] CRAN (R 4.2.2)
# googledrive              2.1.0     2023-03-22 [2] CRAN (R 4.2.3)
# gridExtra                2.3       2017-09-09 [2] CRAN (R 4.2.1)
# gtable                   0.3.3     2023-03-21 [2] CRAN (R 4.2.3)
# here                   * 1.0.1     2020-12-13 [2] CRAN (R 4.2.1)
# hms                      1.1.3     2023-03-21 [2] CRAN (R 4.2.3)
# htmltools                0.5.5     2023-03-23 [2] CRAN (R 4.2.3)
# htmlwidgets              1.6.2     2023-03-17 [2] CRAN (R 4.2.3)
# httpuv                   1.6.9     2023-02-14 [2] CRAN (R 4.2.2)
# httr                     1.4.5     2023-02-24 [2] CRAN (R 4.2.2)
# ica                      1.0-3     2022-07-08 [1] CRAN (R 4.2.2)
# igraph                   1.4.2     2023-04-07 [2] CRAN (R 4.2.3)
# interactiveDisplayBase   1.36.0    2022-11-01 [2] Bioconductor
# IRanges                * 2.32.0    2022-11-01 [2] Bioconductor
# irlba                    2.3.5.1   2022-10-03 [2] CRAN (R 4.2.1)
# jaffelab               * 0.99.32   2023-02-09 [1] Github (LieberInstitute/jaffelab@7b7afe3)
# jsonlite                 1.8.5     2023-06-05 [1] CRAN (R 4.2.3)
# KEGGREST                 1.38.0    2022-11-01 [2] Bioconductor
# KernSmooth               2.23-20   2021-05-03 [3] CRAN (R 4.2.3)
# later                    1.3.0     2021-08-18 [2] CRAN (R 4.2.1)
# lattice                  0.20-45   2021-09-22 [3] CRAN (R 4.2.3)
# lazyeval                 0.2.2     2019-03-15 [2] CRAN (R 4.2.1)
# leiden                   0.4.3     2022-09-10 [1] CRAN (R 4.2.2)
# lifecycle                1.0.3     2022-10-07 [2] CRAN (R 4.2.1)
# limma                    3.54.2    2023-02-28 [2] Bioconductor
# listenv                  0.9.0     2022-12-16 [2] CRAN (R 4.2.2)
# lmtest                   0.9-40    2022-03-21 [2] CRAN (R 4.2.1)
# magrittr                 2.0.3     2022-03-30 [2] CRAN (R 4.2.1)
# MASS                     7.3-58.2  2023-01-23 [3] CRAN (R 4.2.3)
# Matrix                   1.5-4     2023-04-04 [2] CRAN (R 4.2.3)
# MatrixGenerics         * 1.10.0    2022-11-01 [2] Bioconductor
# matrixStats            * 1.0.0     2023-06-02 [1] CRAN (R 4.2.3)
# memoise                  2.0.1     2021-11-26 [2] CRAN (R 4.2.1)
# mime                     0.12      2021-09-28 [2] CRAN (R 4.2.1)
# miniUI                   0.1.1.1   2018-05-18 [2] CRAN (R 4.2.1)
# munsell                  0.5.0     2018-06-12 [2] CRAN (R 4.2.1)
# nlme                     3.1-162   2023-01-31 [2] CRAN (R 4.2.2)
# org.Hs.eg.db           * 3.16.0    2022-09-28 [2] Bioconductor
# parallelly               1.36.0    2023-05-26 [1] CRAN (R 4.2.3)
# patchwork                1.1.2     2022-08-19 [2] CRAN (R 4.2.1)
# pbapply                  1.7-0     2023-01-13 [2] CRAN (R 4.2.2)
# pillar                   1.9.0     2023-03-22 [2] CRAN (R 4.2.3)
# pkgconfig                2.0.3     2019-09-22 [2] CRAN (R 4.2.1)
# plotly                   4.10.1    2022-11-07 [2] CRAN (R 4.2.2)
# plyr                     1.8.8     2022-11-11 [2] CRAN (R 4.2.2)
# png                      0.1-8     2022-11-29 [2] CRAN (R 4.2.2)
# polyclip                 1.10-4    2022-10-20 [2] CRAN (R 4.2.2)
# prettyunits              1.1.1     2020-01-24 [2] CRAN (R 4.2.1)
# progress                 1.2.2     2019-05-16 [2] CRAN (R 4.2.1)
# progressr                0.13.0    2023-01-10 [1] CRAN (R 4.2.2)
# promises                 1.2.0.1   2021-02-11 [2] CRAN (R 4.2.1)
# ProtGenerics             1.30.0    2022-11-01 [2] Bioconductor
# purrr                    1.0.1     2023-01-10 [2] CRAN (R 4.2.2)
# R6                       2.5.1     2021-08-19 [2] CRAN (R 4.2.1)
# rafalib                * 1.0.0     2015-08-09 [1] CRAN (R 4.2.2)
# RANN                     2.6.1     2019-01-08 [2] CRAN (R 4.2.1)
# rappdirs                 0.3.3     2021-01-31 [2] CRAN (R 4.2.1)
# RColorBrewer             1.1-3     2022-04-03 [2] CRAN (R 4.2.1)
# Rcpp                     1.0.10    2023-01-22 [2] CRAN (R 4.2.2)
# RcppAnnoy                0.0.20    2022-10-27 [2] CRAN (R 4.2.2)
# RCurl                    1.98-1.12 2023-03-27 [2] CRAN (R 4.2.3)
# reshape2                 1.4.4     2020-04-09 [2] CRAN (R 4.2.1)
# restfulr                 0.0.15    2022-06-16 [2] CRAN (R 4.2.1)
# reticulate               1.28      2023-01-27 [2] CRAN (R 4.2.2)
# rJava                    1.0-6     2021-12-10 [2] CRAN (R 4.2.1)
# rjson                    0.2.21    2022-01-09 [2] CRAN (R 4.2.1)
# rlang                    1.1.1     2023-04-28 [1] CRAN (R 4.2.3)
# ROCR                     1.0-11    2020-05-02 [2] CRAN (R 4.2.1)
# rprojroot                2.0.3     2022-04-02 [2] CRAN (R 4.2.1)
# Rsamtools                2.14.0    2022-11-01 [2] Bioconductor
# RSQLite                  2.3.1     2023-04-03 [2] CRAN (R 4.2.3)
# rsvd                     1.0.5     2021-04-16 [2] CRAN (R 4.2.1)
# rtracklayer              1.58.0    2022-11-01 [2] Bioconductor
# Rtsne                    0.16      2022-04-17 [2] CRAN (R 4.2.1)
# S4Vectors              * 0.36.2    2023-02-26 [2] Bioconductor
# ScaledMatrix             1.6.0     2022-11-01 [2] Bioconductor
# scales                   1.2.1     2022-08-20 [2] CRAN (R 4.2.1)
# scater                 * 1.26.1    2022-11-13 [2] Bioconductor
# scattermore              1.1       2023-05-17 [1] CRAN (R 4.2.3)
# sctransform              0.3.5     2022-09-21 [1] CRAN (R 4.2.2)
# scuttle                * 1.8.4     2023-01-19 [2] Bioconductor
# segmented                1.6-4     2023-04-13 [1] CRAN (R 4.2.3)
# sessioninfo            * 1.2.2     2021-12-06 [2] CRAN (R 4.2.1)
# Seurat                 * 4.3.0     2022-11-18 [1] CRAN (R 4.2.3)
# SeuratObject           * 4.1.3     2022-11-07 [1] CRAN (R 4.2.2)
# shiny                    1.7.4     2022-12-15 [2] CRAN (R 4.2.2)
# SingleCellExperiment   * 1.20.1    2023-03-17 [2] Bioconductor
# sp                       1.6-1     2023-05-31 [1] CRAN (R 4.2.3)
# sparseMatrixStats        1.10.0    2022-11-01 [2] Bioconductor
# spatstat.data            3.0-1     2023-03-12 [1] CRAN (R 4.2.3)
# spatstat.explore         3.2-1     2023-05-13 [1] CRAN (R 4.2.3)
# spatstat.geom            3.2-1     2023-05-09 [1] CRAN (R 4.2.3)
# spatstat.random          3.1-5     2023-05-11 [1] CRAN (R 4.2.3)
# spatstat.sparse          3.0-1     2023-03-12 [1] CRAN (R 4.2.3)
# spatstat.utils           3.0-3     2023-05-09 [1] CRAN (R 4.2.3)
# stringi                  1.7.12    2023-01-11 [2] CRAN (R 4.2.2)
# stringr                  1.5.0     2022-12-02 [2] CRAN (R 4.2.2)
# SummarizedExperiment   * 1.28.0    2022-11-01 [2] Bioconductor
# survival                 3.5-3     2023-02-12 [3] CRAN (R 4.2.3)
# tensor                   1.5       2012-05-05 [1] CRAN (R 4.2.2)
# tibble                   3.2.1     2023-03-20 [2] CRAN (R 4.2.3)
# tidyr                    1.3.0     2023-01-24 [2] CRAN (R 4.2.2)
# tidyselect               1.2.0     2022-10-10 [2] CRAN (R 4.2.1)
# utf8                     1.2.3     2023-01-31 [2] CRAN (R 4.2.2)
# uwot                     0.1.14    2022-08-22 [2] CRAN (R 4.2.1)
# vctrs                    0.6.2     2023-04-19 [1] CRAN (R 4.2.3)
# vipor                    0.4.5     2017-03-22 [2] CRAN (R 4.2.1)
# viridis                  0.6.2     2021-10-13 [2] CRAN (R 4.2.1)
# viridisLite              0.4.2     2023-05-02 [1] CRAN (R 4.2.3)
# withr                    2.5.0     2022-03-03 [2] CRAN (R 4.2.1)
# xlsx                   * 0.6.5     2020-11-10 [2] CRAN (R 4.2.1)
# xlsxjars                 0.6.1     2014-08-22 [2] CRAN (R 4.2.1)
# XML                      3.99-0.14 2023-03-19 [2] CRAN (R 4.2.3)
# xml2                     1.3.3     2021-11-30 [2] CRAN (R 4.2.1)
# xtable                   1.8-4     2019-04-21 [2] CRAN (R 4.2.1)
# XVector                  0.38.0    2022-11-01 [2] Bioconductor
# yaml                     2.3.7     2023-01-23 [2] CRAN (R 4.2.2)
# zlibbioc                 1.44.0    2022-11-01 [2] Bioconductor
# zoo                      1.8-12    2023-04-13 [1] CRAN (R 4.2.3)
# 
# [1] /users/bsimbiat/R/4.2.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/library
