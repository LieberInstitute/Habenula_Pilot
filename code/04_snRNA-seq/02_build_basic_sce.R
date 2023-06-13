# January 6, 2023 - Bukola Ajanaku
# Building basic sce objects. This follows the initial pca analysis for bulk 
# Habenula data. Only 7 samples, all control.
# Based on: 
# 1) https://github.com/LieberInstitute/Habenula_Bulk/blob/master/code/09_snRNA-seq_re-processed/02_qc.R
# 2) https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/code/03_build_sce/01_build_basic_sce.R
# 3) https://github.com/LieberInstitute/Habenula_Bulk/blob/master/code/09_snRNA-seq_re-processed/00_get_droplet_scores.R
# qrsh -l mem_free=75G,h_vmem=75G

# Loading relevant libraries
library("SingleCellExperiment")
library("SpatialExperiment")
library("DropletUtils")
library("here")
library("rtracklayer")
library("lobstr")
library("sessioninfo")
library("dplyr")
library("scater")
library("purrr")

# Sample names 
sample_list <- list.files(here("processed-data", "07_cellranger"))

sce_list <- map(sample_list, function(sample){
  fn10x <-  here("processed-data", "07_cellranger", sample, "outs", "raw_feature_bc_matrix") 
  fnDropScore <- here("processed-data", "04_snRNA-seq", "01_get_droplet_scores", paste0("droplet_scores_", sample, ".Rdata"))
  
  sce <- read10xCounts(fn10x, col.names=TRUE)
  load(fnDropScore, verbose = TRUE)
  ncol_preDrop = ncol(sce)
  # 1148322
  
  sce <- sce[, which(e.out$FDR <= 0.001)]
  ncol_postDrop = ncol(sce)
  # 3622
  
  message(sample, ": Pre-Drop = ", ncol_preDrop, " and Post-Drop = ", ncol_postDrop)

  return(sce)
})

# Match rownames and appending by column
make_Rows <- rownames(sce_list[[1]])
sce_list <- map(sce_list, ~.x[make_Rows,])
# identical(rownames(sce_list2[[1]]), rownames(sce_list2[[7]]))
# TRUE

# OFFICIAL SCE HB OBJECT PRE-QC. 
sce_hb_preQC <- do.call("cbind", sce_list)

save(sce_hb_preQC, file = here("processed-data", "04_snRNA-seq", "sce_objects", 
                               "sce_hb_preQC.Rdata"))
# sgejobs::job_single('02_build_basic_sce', create_shell = TRUE, queue= 'bluejay', memory = '100G', command = "Rscript 02_build_basic_sce.R")


## QUALITY CONTROL
# Thought process: (this is 10x Genomics data)
# 1) Based on OSCA book, can use mean absolute deviation (MADs) approach to 
# determine outliers but this is not as "straightforward" as it may through 
# away neurons.
# 2) There is a test using UMI/barcode rank (the knee plots), to determine 
# individual thresholds for what is too low of quality but this may drop cells 
# with organically low RNA content.
# 3) This is a single-nuc RNA project meaning that dropping high mito content
# is a useful QC metric because it shows samples where cytoplasm 
# was not fully or successfully stripped.
# Game plan:
# 1) Drop empty droplets. (script 1, completed)
# 2) Build sce objects. (this script, script 2)
# 3) Dropping high mito content (script 3, next script)


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
# beachmat               2.14.2    2023-04-07 [2] Bioconductor
# beeswarm               0.4.0     2021-06-01 [2] CRAN (R 4.2.1)
# Biobase              * 2.58.0    2022-11-01 [2] Bioconductor
# BiocGenerics         * 0.44.0    2022-11-01 [2] Bioconductor
# BiocIO                 1.8.0     2022-11-01 [2] Bioconductor
# BiocNeighbors          1.16.0    2022-11-01 [2] Bioconductor
# BiocParallel           1.32.6    2023-03-17 [2] Bioconductor
# BiocSingular           1.14.0    2022-11-01 [2] Bioconductor
# Biostrings             2.66.0    2022-11-01 [2] Bioconductor
# bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.2.1)
# cli                    3.6.1     2023-03-23 [2] CRAN (R 4.2.3)
# codetools              0.2-19    2023-02-01 [3] CRAN (R 4.2.3)
# colorout             * 1.2-2     2023-02-09 [1] Github (jalvesaq/colorout@79931fd)
# colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.2.2)
# crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.2.1)
# DelayedArray           0.24.0    2022-11-01 [2] Bioconductor
# DelayedMatrixStats     1.20.0    2022-11-01 [2] Bioconductor
# dplyr                * 1.1.1     2023-03-22 [2] CRAN (R 4.2.3)
# dqrng                  0.3.0     2021-05-01 [2] CRAN (R 4.2.1)
# DropletUtils         * 1.18.1    2022-11-22 [2] Bioconductor
# edgeR                  3.40.2    2023-01-19 [2] Bioconductor
# fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.2.2)
# generics               0.1.3     2022-07-05 [2] CRAN (R 4.2.1)
# GenomeInfoDb         * 1.34.9    2023-02-02 [2] Bioconductor
# GenomeInfoDbData       1.2.9     2022-09-29 [2] Bioconductor
# GenomicAlignments      1.34.1    2023-03-09 [2] Bioconductor
# GenomicRanges        * 1.50.2    2022-12-16 [2] Bioconductor
# ggbeeswarm             0.7.1     2022-12-16 [2] CRAN (R 4.2.2)
# ggplot2              * 3.4.2     2023-04-03 [2] CRAN (R 4.2.3)
# ggrepel                0.9.3     2023-02-03 [2] CRAN (R 4.2.2)
# glue                   1.6.2     2022-02-24 [2] CRAN (R 4.2.1)
# gridExtra              2.3       2017-09-09 [2] CRAN (R 4.2.1)
# gtable                 0.3.3     2023-03-21 [2] CRAN (R 4.2.3)
# HDF5Array              1.26.0    2022-11-01 [2] Bioconductor
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.2.1)
# IRanges              * 2.32.0    2022-11-01 [2] Bioconductor
# irlba                  2.3.5.1   2022-10-03 [2] CRAN (R 4.2.1)
# lattice                0.20-45   2021-09-22 [3] CRAN (R 4.2.3)
# lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.2.1)
# limma                  3.54.2    2023-02-28 [2] Bioconductor
# lobstr               * 1.1.2     2022-06-22 [2] CRAN (R 4.2.1)
# locfit                 1.5-9.7   2023-01-02 [2] CRAN (R 4.2.2)
# magick                 2.7.4     2023-03-09 [2] CRAN (R 4.2.3)
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.2.1)
# Matrix                 1.5-4     2023-04-04 [2] CRAN (R 4.2.3)
# MatrixGenerics       * 1.10.0    2022-11-01 [2] Bioconductor
# matrixStats          * 1.0.0     2023-06-02 [1] CRAN (R 4.2.3)
# munsell                0.5.0     2018-06-12 [2] CRAN (R 4.2.1)
# pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.2.3)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.2.1)
# purrr                * 1.0.1     2023-01-10 [2] CRAN (R 4.2.2)
# R.methodsS3            1.8.2     2022-06-13 [2] CRAN (R 4.2.1)
# R.oo                   1.25.0    2022-06-12 [2] CRAN (R 4.2.1)
# R.utils                2.12.2    2022-11-11 [2] CRAN (R 4.2.2)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.2.1)
# Rcpp                   1.0.10    2023-01-22 [2] CRAN (R 4.2.2)
# RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.2.3)
# restfulr               0.0.15    2022-06-16 [2] CRAN (R 4.2.1)
# rhdf5                  2.42.1    2023-04-07 [2] Bioconductor
# rhdf5filters           1.10.1    2023-03-24 [2] Bioconductor
# Rhdf5lib               1.20.0    2022-11-01 [2] Bioconductor
# rjson                  0.2.21    2022-01-09 [2] CRAN (R 4.2.1)
# rlang                  1.1.1     2023-04-28 [1] CRAN (R 4.2.3)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.2.1)
# Rsamtools              2.14.0    2022-11-01 [2] Bioconductor
# rsvd                   1.0.5     2021-04-16 [2] CRAN (R 4.2.1)
# rtracklayer          * 1.58.0    2022-11-01 [2] Bioconductor
# S4Vectors            * 0.36.2    2023-02-26 [2] Bioconductor
# ScaledMatrix           1.6.0     2022-11-01 [2] Bioconductor
# scales                 1.2.1     2022-08-20 [2] CRAN (R 4.2.1)
# scater               * 1.26.1    2022-11-13 [2] Bioconductor
# scuttle              * 1.8.4     2023-01-19 [2] Bioconductor
# sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.2.1)
# SingleCellExperiment * 1.20.1    2023-03-17 [2] Bioconductor
# sparseMatrixStats      1.10.0    2022-11-01 [2] Bioconductor
# SpatialExperiment    * 1.8.1     2023-03-05 [2] Bioconductor
# SummarizedExperiment * 1.28.0    2022-11-01 [2] Bioconductor
# tibble                 3.2.1     2023-03-20 [2] CRAN (R 4.2.3)
# tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.2.1)
# utf8                   1.2.3     2023-01-31 [2] CRAN (R 4.2.2)
# vctrs                  0.6.2     2023-04-19 [1] CRAN (R 4.2.3)
# vipor                  0.4.5     2017-03-22 [2] CRAN (R 4.2.1)
# viridis                0.6.2     2021-10-13 [2] CRAN (R 4.2.1)
# viridisLite            0.4.2     2023-05

