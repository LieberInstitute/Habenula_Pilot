## May 2, 2023 - Bukola Ajanaku
# sce data pre-qc on violin plots to show threshold for drops. 
# qrsh -l mem_free=30G,h_vmem=30G

# loading relevant libraries
library("SingleCellExperiment")
library("here")
library("ggplot2")
library("EnsDb.Hsapiens.v86")
library("scuttle")
library("jaffelab")
library("sessioninfo")
library("scater")
library("PupillometryR")
library("cowplot")

# loading pre-QC sce object
load(here("processed-data", "99_paper_figs",  "sce_objects", "pre_QC_sce.Rdata"),
     verbose = TRUE)
# sce_hb_preQC
sce <- sce_hb_preQC
rm(sce_hb_preQC)

# creating plot directory
plot_dir <- here("plots", "99_paper_figs", "04_QC_Violin.R")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

############### FINDING HIGH MITO ##############################################
# MADs approach for high mito droplets (indicates nuclei where cytoplasm wasn't
# successfully stripped):
set.seed(777)

# Defining locations
location <- mapIds(EnsDb.Hsapiens.v86, keys = rowData(sce)$ID, 
                   column = "SEQNAME", keytype = "GENEID")
# Unable to map 2611 of 36601 requested IDs.

# Identifying mito reads
stats <- perCellQCMetrics(sce, subsets = list(Mito = which(location=="MT")))

# Changing Sample ID names from locations 
sce$path <- sce$Sample
sce$Sample <- ss(sce$Sample, "/", 9)

# Binding stats to colData of sce
colData(sce) <- cbind(colData(sce), stats)

############### HIGH MITO ######################################################
sce$high_mito <- isOutlier(sce$subsets_Mito_percent, nmads=3, type="higher", 
                           batch = sce$Sample)
table(sce$high_mito, by = sce$Sample)
#       Br1092 Br1204 Br1469 Br1735 Br5555 Br5558 Br5639
# FALSE   3352   1391   2185   3199   3600    801   3174
# TRUE     270    274    314    544    305    135    258

############### LOW LIBRARY SIZE ###############################################
sce$lowLib <- isOutlier(sce$sum, log=TRUE, type="lower", batch = sce$Sample)
table(sce$lowLib, by = sce$Sample)
#       Br1092 Br1204 Br1469 Br1735 Br5555 Br5558 Br5639
# FALSE   3430   1445   2414   3633   3772    847   3348
# TRUE     192    220     85    110    133     89     84

############### LOW DETECTED FEATURES ##########################################
sce$lowDetecFea <- isOutlier(sce$detected, log=TRUE, type="lower", batch = sce$Sample)
table(sce$lowDetecFea, by = sce$Sample)
#       Br1092 Br1204 Br1469 Br1735 Br5555 Br5558 Br5639
# FALSE   3394   1334   2347   3437   3766    826   3289
# TRUE     228    331    152    306    139    110    143


################ PLOTTING PER METRIC ###########################################
pd <- as.data.frame(colData(sce))
sce$lowLib <- as.factor(sce$lowLib)
sce$Sample <- as.factor(sce$Sample)

## high_mito
a <-  plotColData(sce, x = "Sample", y="subsets_Mito_percent", colour_by="high_mito") +
  scale_y_log10() + 
  labs(y = "Mito Percent") + 
  scale_colour_discrete(name="High Mito?") +
  aes(group = pd$Sample) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_blank())

pdf(here(plot_dir, "official_Violin_QC_High_Mito.pdf"), width = 8)
  a
dev.off()

# low libray size
b <-  plotColData(sce, x = "Sample", y="sum", colour_by="lowLib") +
  scale_y_log10() + 
  labs(y = "Library Size") + 
  scale_colour_discrete(name="Low Lib Size?") +
  aes(group = pd$Sample) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_blank())

pdf(here(plot_dir, "official_Violin_QC_Low_Lib.pdf"), width = 8)
  b
dev.off()


# low detected features
c <-  plotColData(sce, x = "Sample", y="detected", colour_by="lowDetecFea") +
  scale_y_log10() + 
  labs(y = "Detected Feature") + 
  scale_colour_discrete(name="Low Detected?") +
  aes(group = pd$Sample) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_blank())

pdf(here(plot_dir, "official_Violin_QC_Low_Det_Feat.pdf"), width = 8)
  c
dev.off()

pdf(here(plot_dir, "official_All_QC_Violin_Plots.pdf"), width = 10, height = 13)
  a <- a + theme(axis.text.x = element_blank())
  b <- b + theme(axis.text.x = element_blank())
  plot_grid(a, b, c, ncol = 1)
dev.off()

####### FOR ONE DRIVE ##########################################################
## high_mito

pdf(here(plot_dir, "forOneDrive", "sfigu_mito_percent_qc_violin.pdf"), 
    height = 5, width = 7.5)

      plotColData(sce, x = "Sample", y="subsets_Mito_percent", colour_by="high_mito") +
        scale_y_log10() + 
        labs(y = "Mito Percent") + 
        scale_colour_discrete(name="High Mito?") +
        aes(group = pd$Sample) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), 
              axis.title.x = element_blank())
dev.off()

# low libray size
pdf(here(plot_dir, "forOneDrive", "sfigu_lib_size_qc_violin.pdf"), 
    height = 5, width = 7.5)

      plotColData(sce, x = "Sample", y="sum", colour_by="lowLib") +
        scale_y_log10() + 
        labs(y = "Library Size") + 
        scale_colour_discrete(name="Low Lib Size?") +
        aes(group = pd$Sample) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), 
              axis.title.x = element_blank())
dev.off()


# low detected features
pdf(here(plot_dir, "forOneDrive" , "sfigu_detected_feats_qc_violin.pdf"), 
      height = 5, width = 7.5)

        plotColData(sce, x = "Sample", y="detected", colour_by="lowDetecFea") +
          scale_y_log10() + 
          labs(y = "Detected Features") + 
          scale_colour_discrete(name="Low Detected?") +
          aes(group = pd$Sample) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1), 
                axis.title.x = element_blank())
dev.off()

sessioninfo::session_info()

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
# package              * version   date (UTC) lib source
# AnnotationDbi        * 1.60.2    2023-03-10 [2] Bioconductor
# AnnotationFilter     * 1.22.0    2022-11-01 [2] Bioconductor
# beachmat               2.14.2    2023-04-07 [2] Bioconductor
# beeswarm               0.4.0     2021-06-01 [2] CRAN (R 4.2.1)
# Biobase              * 2.58.0    2022-11-01 [2] Bioconductor
# BiocFileCache          2.6.1     2023-02-17 [2] Bioconductor
# BiocGenerics         * 0.44.0    2022-11-01 [2] Bioconductor
# BiocIO                 1.8.0     2022-11-01 [2] Bioconductor
# BiocNeighbors          1.16.0    2022-11-01 [2] Bioconductor
# BiocParallel           1.32.6    2023-03-17 [2] Bioconductor
# BiocSingular           1.14.0    2022-11-01 [2] Bioconductor
# biomaRt                2.54.1    2023-03-20 [2] Bioconductor
# Biostrings             2.66.0    2022-11-01 [2] Bioconductor
# bit                    4.0.5     2022-11-15 [2] CRAN (R 4.2.2)
# bit64                  4.0.5     2020-08-30 [2] CRAN (R 4.2.1)
# bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.2.1)
# blob                   1.2.4     2023-03-17 [2] CRAN (R 4.2.3)
# cachem                 1.0.8     2023-05-01 [1] CRAN (R 4.2.3)
# cli                    3.6.1     2023-03-23 [2] CRAN (R 4.2.3)
# codetools              0.2-19    2023-02-01 [3] CRAN (R 4.2.3)
# colorout             * 1.2-2     2023-02-09 [1] Github (jalvesaq/colorout@79931fd)
# colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.2.2)
# cowplot              * 1.1.1     2020-12-30 [2] CRAN (R 4.2.1)
# crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.2.1)
# curl                   5.0.0     2023-01-12 [2] CRAN (R 4.2.2)
# DBI                    1.1.3     2022-06-18 [2] CRAN (R 4.2.1)
# dbplyr                 2.3.2     2023-03-21 [2] CRAN (R 4.2.3)
# DelayedArray           0.24.0    2022-11-01 [2] Bioconductor
# DelayedMatrixStats     1.20.0    2022-11-01 [2] Bioconductor
# digest                 0.6.31    2022-12-11 [2] CRAN (R 4.2.2)
# dplyr                * 1.1.1     2023-03-22 [2] CRAN (R 4.2.3)
# EnsDb.Hsapiens.v86   * 2.99.0    2023-02-09 [1] Bioconductor
# ensembldb            * 2.22.0    2022-11-01 [2] Bioconductor
# fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.2.2)
# fastmap                1.1.1     2023-02-24 [2] CRAN (R 4.2.2)
# filelock               1.0.2     2018-10-05 [2] CRAN (R 4.2.1)
# fs                     1.6.1     2023-02-06 [2] CRAN (R 4.2.2)
# gargle                 1.3.0     2023-01-30 [2] CRAN (R 4.2.2)
# generics               0.1.3     2022-07-05 [2] CRAN (R 4.2.1)
# GenomeInfoDb         * 1.34.9    2023-02-02 [2] Bioconductor
# GenomeInfoDbData       1.2.9     2022-09-29 [2] Bioconductor
# GenomicAlignments      1.34.1    2023-03-09 [2] Bioconductor
# GenomicFeatures      * 1.50.4    2023-01-24 [2] Bioconductor
# GenomicRanges        * 1.50.2    2022-12-16 [2] Bioconductor
# ggbeeswarm             0.7.1     2022-12-16 [2] CRAN (R 4.2.2)
# ggplot2              * 3.4.2     2023-04-03 [2] CRAN (R 4.2.3)
# ggrepel                0.9.3     2023-02-03 [2] CRAN (R 4.2.2)
# glue                   1.6.2     2022-02-24 [2] CRAN (R 4.2.1)
# googledrive            2.1.0     2023-03-22 [2] CRAN (R 4.2.3)
# gridExtra              2.3       2017-09-09 [2] CRAN (R 4.2.1)
# gtable                 0.3.3     2023-03-21 [2] CRAN (R 4.2.3)
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.2.1)
# hms                    1.1.3     2023-03-21 [2] CRAN (R 4.2.3)
# httr                   1.4.5     2023-02-24 [2] CRAN (R 4.2.2)
# IRanges              * 2.32.0    2022-11-01 [2] Bioconductor
# irlba                  2.3.5.1   2022-10-03 [2] CRAN (R 4.2.1)
# jaffelab             * 0.99.32   2023-02-09 [1] Github (LieberInstitute/jaffelab@7b7afe3)
# KEGGREST               1.38.0    2022-11-01 [2] Bioconductor
# lattice                0.20-45   2021-09-22 [3] CRAN (R 4.2.3)
# lazyeval               0.2.2     2019-03-15 [2] CRAN (R 4.2.1)
# lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.2.1)
# limma                  3.54.2    2023-02-28 [2] Bioconductor
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.2.1)
# MASS                   7.3-58.2  2023-01-23 [3] CRAN (R 4.2.3)
# Matrix                 1.5-4     2023-04-04 [2] CRAN (R 4.2.3)
# MatrixGenerics       * 1.10.0    2022-11-01 [2] Bioconductor
# matrixStats          * 1.0.0     2023-06-02 [1] CRAN (R 4.2.3)
# memoise                2.0.1     2021-11-26 [2] CRAN (R 4.2.1)
# munsell                0.5.0     2018-06-12 [2] CRAN (R 4.2.1)
# nlme                   3.1-162   2023-01-31 [2] CRAN (R 4.2.2)
# pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.2.3)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.2.1)
# png                    0.1-8     2022-11-29 [2] CRAN (R 4.2.2)
# prettyunits            1.1.1     2020-01-24 [2] CRAN (R 4.2.1)
# progress               1.2.2     2019-05-16 [2] CRAN (R 4.2.1)
# ProtGenerics           1.30.0    2022-11-01 [2] Bioconductor
# PupillometryR        * 0.0.4     2021-09-19 [1] CRAN (R 4.2.3)
# purrr                  1.0.1     2023-01-10 [2] CRAN (R 4.2.2)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.2.1)
# rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 4.2.2)
# rappdirs               0.3.3     2021-01-31 [2] CRAN (R 4.2.1)
# RColorBrewer           1.1-3     2022-04-03 [2] CRAN (R 4.2.1)
# Rcpp                   1.0.10    2023-01-22 [2] CRAN (R 4.2.2)
# RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.2.3)
# restfulr               0.0.15    2022-06-16 [2] CRAN (R 4.2.1)
# rjson                  0.2.21    2022-01-09 [2] CRAN (R 4.2.1)
# rlang                * 1.1.1     2023-04-28 [1] CRAN (R 4.2.3)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.2.1)
# Rsamtools              2.14.0    2022-11-01 [2] Bioconductor
# RSQLite                2.3.1     2023-04-03 [2] CRAN (R 4.2.3)
# rsvd                   1.0.5     2021-04-16 [2] CRAN (R 4.2.1)
# rtracklayer            1.58.0    2022-11-01 [2] Bioconductor
# S4Vectors            * 0.36.2    2023-02-26 [2] Bioconductor
# ScaledMatrix           1.6.0     2022-11-01 [2] Bioconductor
# scales                 1.2.1     2022-08-20 [2] CRAN (R 4.2.1)
# scater               * 1.26.1    2022-11-13 [2] Bioconductor
# scuttle              * 1.8.4     2023-01-19 [2] Bioconductor
# segmented              1.6-4     2023-04-13 [1] CRAN (R 4.2.3)
# sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.2.1)
# SingleCellExperiment * 1.20.1    2023-03-17 [2] Bioconductor
# sparseMatrixStats      1.10.0    2022-11-01 [2] Bioconductor
# stringi                1.7.12    2023-01-11 [2] CRAN (R 4.2.2)
# stringr                1.5.0     2022-12-02 [2] CRAN (R 4.2.2)
# SummarizedExperiment * 1.28.0    2022-11-01 [2] Bioconductor
# tibble                 3.2.1     2023-03-20 [2] CRAN (R 4.2.3)
# tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.2.1)
# utf8                   1.2.3     2023-01-31 [2] CRAN (R 4.2.2)
# vctrs                  0.6.2     2023-04-19 [1] CRAN (R 4.2.3)
# vipor                  0.4.5     2017-03-22 [2] CRAN (R 4.2.1)
# viridis                0.6.2     2021-10-13 [2] CRAN (R 4.2.1)
# viridisLite            0.4.2     2023-05-02 [1] CRAN (R 4.2.3)
# withr                  2.5.0     2022-03-03 [2] CRAN (R 4.2.1)
# XML                    3.99-0.14 2023-03-19
