# January 18, 2023 - Bukola Ajanaku
# Normalizing sce object by different metrics and then plotting PCs, TSNEs, and
# UMAPs.
# Based on:
# 1) https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/code/03_build_sce/05_harmony_correction.R
# 2) https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/code/03_build_sce/05.5_harmony_correction_plots.R
# qrsh -l mem_free=50G,h_vmem=50G

# Loading relevant libraries
library("SingleCellExperiment")
library("harmony")
library("scater")
library("here")
library("sessioninfo")
library("HDF5Array")
library("viridis")
library("ggplot2")

# Loading sce_uncorrected 
# if locations ever stop working, please check 99_paper_figs folder in code 
# to search for new possible location
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
          "sce_uncorrected_PCA.Rdata"))
sce_uncorrected

# BUG FIX: Uniquifying colnames using barcode and sample #######################
# BUG:
table(duplicated(colnames(sce_uncorrected)))
# FALSE  TRUE 
# 17048    34 

# test <- "TCGCTCACAAATAGCA-1"
# colData(sce)[sce$Barcode == test,]
# DataFrame with 2 rows and 16 columns
# Sample            Barcode                   path
# <character>        <character>            <character> 
#   TCGCTCACAAATAGCA-1      Br1092 TCGCTCACAAATAGCA-1 /dcs04/lieber/lcolla..
# TCGCTCACAAATAGCA-1      Br1204 TCGCTCACAAATAGCA-1 /dcs04/lieber/lcolla..

colnames(sce_uncorrected) <- paste0(sce_uncorrected$Sample, "_", 
                                    sce_uncorrected$Barcode)
any(duplicated(colnames(sce_uncorrected)))
# [1] FALSE 


####### RUNNING HARMONY ########################################################
# Note: We used GLM PCA not normal PCA but harmony function searches for something
# named "PCA". Thus, we must submit GLM under the name PCA. After running harmony,
# we can do away with our fake PCA:

# Submitting GLM PCA under the name PCA:
reducedDim(sce_uncorrected, "PCA") <- reducedDim(sce_uncorrected, "GLMPCA_approx")

# Running harmony
sce_corrbySamp <- RunHarmony(sce_uncorrected, group.by.vars = "Sample", verbose = TRUE)
# sce_corrbyRun <- RunHarmony(sce_uncorrected, group.by.vars = "Run", verbose = TRUE)
    # Error in harmonyObj$init_cluster_cpp(0) : 
    # element-wise multiplication: incompatible matrix dimensions: 100x3 and 100x1

# Removing our redudant fake PCA name and keeping it under "GLMPCA_approx"
reducedDim(sce_uncorrected, "PCA") <- NULL

####### TSNE & UMAP ############################################################
set.seed(777)

## sce corrected by Sample
sce_corrbySamp <- runTSNE(sce_corrbySamp, dimred = "HARMONY")
sce_corrbySamp <- runUMAP(sce_corrbySamp, dimred = "HARMONY")

####### PLOTTING ###############################################################
# GLM by Sample
pdf(here("plots", "04_snRNA-seq", "05_GLM_Harmony_plots", "GLM_harmony_plot_by_Sample.pdf"))
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "Sample") 
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "Sample") + 
    facet_wrap(~ sce_corrbySamp$Sample)
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "Sample", ncomponents = c(2,3)) # looks batchy
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "Sample", ncomponents = c(2,3)) + 
    facet_wrap(~ sce_corrbySamp$Sample)
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "Sample", ncomponents = 5)
dev.off()

# Plotting using GLM-PCA functions [color by Run]
pdf(here("plots", "04_snRNA-seq", "05_GLM_Harmony_plots", "GLM_harmony_plot_by_Run.pdf"))
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "Run") 
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "Run") + 
    facet_wrap(~ sce_corrbySamp$Sample)              
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "Run", ncomponents = c(2,3)) 
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "Run", ncomponents = c(2,3)) + 
    facet_wrap(~ sce_corrbySamp$Sample)
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "Run", ncomponents = 5)
dev.off()

# Plotting using GLM-PCA functions [color by Erik Cluster]
pdf(here("plots", "04_snRNA-seq", "05_GLM_Harmony_plots", "GLM_harmony_plot_by_ct_Erik.pdf"))
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "ct_Erik") + labs(caption = "887 NAs")
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "ct_Erik") + 
    facet_wrap(~ sce_corrbySamp$Sample)  + labs(caption = "887 NAs")            
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "ct_Erik", ncomponents = c(2,3)) + 
    labs(caption = "887 NAs")
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "ct_Erik", ncomponents = c(2,3)) + 
    facet_wrap(~ sce_corrbySamp$Sample) + labs(caption = "887 NAs")
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "ct_Erik", ncomponents = 5) + 
    labs(caption = "887 NAs")
dev.off()

# GLM by Continous Metrics 
pdf(here("plots", "04_snRNA-seq", "05_GLM_Harmony_plots", "GLM_harmony_continous_metrics.pdf"))
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "sum") +
    scale_color_viridis(option = "G", begin = 0.19) 
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "sum")  + 
    facet_wrap(~ sce_corrbySamp$Sample) +
    scale_color_viridis(option = "G", begin = 0.19) 
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "detected") +
    scale_color_viridis(option = "G", begin = 0.19) 
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "detected")  + 
    facet_wrap(~ sce_corrbySamp$Sample) +
    scale_color_viridis(option = "G", begin = 0.19) 
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "subsets_Mito_percent") +
    scale_color_viridis(option = "G", begin = 0.19) 
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "subsets_Mito_percent") + 
    facet_wrap(~ sce_corrbySamp$Sample) +
    scale_color_viridis(option = "G", begin = 0.19) 
dev.off()

# TSNE by Sample
pdf(here("plots", "04_snRNA-seq", "05_GLM_Harmony_plots", "TSNE_harmony_plot.pdf"))
  plotReducedDim(sce_corrbySamp, dimred = "TSNE", colour_by = "Sample") 
  plotReducedDim(sce_corrbySamp, dimred = "TSNE", colour_by = "Sample") + facet_wrap(~ sce_corrbySamp$Sample)
  plotReducedDim(sce_corrbySamp, dimred = "TSNE", colour_by = "ct_Erik") 
  plotReducedDim(sce_corrbySamp, dimred = "TSNE", colour_by = "ct_Erik") + facet_wrap(~ sce_corrbySamp$Sample)
  plotReducedDim(sce_corrbySamp, dimred = "TSNE", colour_by = "Run") 
  plotReducedDim(sce_corrbySamp, dimred = "TSNE", colour_by = "Run") + facet_wrap(~ sce_corrbySamp$Sample)
dev.off()

# TSNE by continuous data
pdf(here("plots", "04_snRNA-seq", "05_GLM_Harmony_plots", "TSNE_harmony_plot_continous_mets.pdf"))
  plotReducedDim(sce_corrbySamp, dimred = "TSNE", colour_by = "sum") 
  plotReducedDim(sce_corrbySamp, dimred = "TSNE", colour_by = "sum") + facet_wrap(~ sce_corrbySamp$Sample)
  plotReducedDim(sce_corrbySamp, dimred = "TSNE", colour_by = "detected") 
  plotReducedDim(sce_corrbySamp, dimred = "TSNE", colour_by = "detected") + facet_wrap(~ sce_corrbySamp$Sample)
  plotReducedDim(sce_corrbySamp, dimred = "TSNE", colour_by = "subsets_Mito_percent") 
  plotReducedDim(sce_corrbySamp, dimred = "TSNE", colour_by = "subsets_Mito_percent") + facet_wrap(~ sce_corrbySamp$Sample)
dev.off()

# UMAP by Sample 
pdf(here("plots", "04_snRNA-seq", "05_GLM_Harmony_plots", "UMAP_harmony_plot.pdf"))
  plotReducedDim(sce_corrbySamp, dimred = "UMAP", colour_by = "Sample")
  plotReducedDim(sce_corrbySamp, dimred = "UMAP", colour_by = "Sample") + facet_wrap(~ sce_corrbySamp$Sample)
  plotReducedDim(sce_corrbySamp, dimred = "UMAP", colour_by = "ct_Erik") 
  plotReducedDim(sce_corrbySamp, dimred = "UMAP", colour_by = "ct_Erik") + facet_wrap(~ sce_corrbySamp$Sample)
  plotReducedDim(sce_corrbySamp, dimred = "UMAP", colour_by = "Run") 
  plotReducedDim(sce_corrbySamp, dimred = "UMAP", colour_by = "Run") + facet_wrap(~ sce_corrbySamp$Sample)
dev.off()

# UMAP by continuous data 
pdf(here("plots", "04_snRNA-seq", "05_GLM_Harmony_plots", "UMAP_harmony_plot_continous_mets.pdf"))
  plotReducedDim(sce_corrbySamp, dimred = "UMAP", colour_by = "sum") 
  plotReducedDim(sce_corrbySamp, dimred = "UMAP", colour_by = "sum") + facet_wrap(~ sce_corrbySamp$Sample)
  plotReducedDim(sce_corrbySamp, dimred = "UMAP", colour_by = "detected") 
  plotReducedDim(sce_corrbySamp, dimred = "UMAP", colour_by = "detected") + facet_wrap(~ sce_corrbySamp$Sample)
  plotReducedDim(sce_corrbySamp, dimred = "UMAP", colour_by = "subsets_Mito_percent") 
  plotReducedDim(sce_corrbySamp, dimred = "UMAP", colour_by = "subsets_Mito_percent") + facet_wrap(~ sce_corrbySamp$Sample)
dev.off()

####### saving ###############################################################
## Saving harmonized (by Sample) sce object 
save(sce_corrbySamp, file = here("processed-data", "04_snRNA-seq", "sce_objects", 
                                  "sce_harmony_by_Samp.Rdata"))

# Saving harmonized (by Sample) sce object (Saved as HDF5 for later clustering)
saveHDF5SummarizedExperiment(sce_corrbySamp, dir = here("processed-data", "04_snRNA-seq", 
                                  "sce_objects", "sce_harmony_by_Samp"), replace = TRUE)

# Reproducibility Information
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
# BiocNeighbors          1.16.0    2022-11-01 [2] Bioconductor
# BiocParallel           1.32.6    2023-03-17 [2] Bioconductor
# BiocSingular           1.14.0    2022-11-01 [2] Bioconductor
# bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.2.1)
# cli                    3.6.1     2023-03-23 [2] CRAN (R 4.2.3)
# codetools              0.2-19    2023-02-01 [3] CRAN (R 4.2.3)
# colorout             * 1.2-2     2023-02-09 [1] Github (jalvesaq/colorout@79931fd)
# colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.2.2)
# cowplot                1.1.1     2020-12-30 [2] CRAN (R 4.2.1)
# DelayedArray         * 0.24.0    2022-11-01 [2] Bioconductor
# DelayedMatrixStats     1.20.0    2022-11-01 [2] Bioconductor
# dplyr                  1.1.1     2023-03-22 [2] CRAN (R 4.2.3)
# fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.2.2)
# generics               0.1.3     2022-07-05 [2] CRAN (R 4.2.1)
# GenomeInfoDb         * 1.34.9    2023-02-02 [2] Bioconductor
# GenomeInfoDbData       1.2.9     2022-09-29 [2] Bioconductor
# GenomicRanges        * 1.50.2    2022-12-16 [2] Bioconductor
# ggbeeswarm             0.7.1     2022-12-16 [2] CRAN (R 4.2.2)
# ggplot2              * 3.4.2     2023-04-03 [2] CRAN (R 4.2.3)
# ggrepel                0.9.3     2023-02-03 [2] CRAN (R 4.2.2)
# glue                   1.6.2     2022-02-24 [2] CRAN (R 4.2.1)
# gridExtra              2.3       2017-09-09 [2] CRAN (R 4.2.1)
# gtable                 0.3.3     2023-03-21 [2] CRAN (R 4.2.3)
# harmony              * 0.1.1     2022-11-14 [2] CRAN (R 4.2.2)
# HDF5Array            * 1.26.0    2022-11-01 [2] Bioconductor
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.2.1)
# IRanges              * 2.32.0    2022-11-01 [2] Bioconductor
# irlba                  2.3.5.1   2022-10-03 [2] CRAN (R 4.2.1)
# lattice                0.20-45   2021-09-22 [3] CRAN (R 4.2.3)
# lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.2.1)
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.2.1)
# Matrix               * 1.5-4     2023-04-04 [2] CRAN (R 4.2.3)
# MatrixGenerics       * 1.10.0    2022-11-01 [2] Bioconductor
# matrixStats          * 1.0.0     2023-06-02 [1] CRAN (R 4.2.3)
# munsell                0.5.0     2018-06-12 [2] CRAN (R 4.2.1)
# pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.2.3)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.2.1)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.2.1)
# Rcpp                 * 1.0.10    2023-01-22 [2] CRAN (R 4.2.2)
# RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.2.3)
# rhdf5                * 2.42.1    2023-04-07 [2] Bioconductor
# rhdf5filters           1.10.1    2023-03-24 [2] Bioconductor
# Rhdf5lib               1.20.0    2022-11-01 [2] Bioconductor
# rlang                  1.1.1     2023-04-28 [1] CRAN (R 4.2.3)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.2.1)
# rsvd                   1.0.5     2021-04-16 [2] CRAN (R 4.2.1)
# S4Vectors            * 0.36.2    2023-02-26 [2] Bioconductor
# ScaledMatrix           1.6.0     2022-11-01 [2] Bioconductor
# scales                 1.2.1     2022-08-20 [2] CRAN (R 4.2.1)
# scater               * 1.26.1    2022-11-13 [2] Bioconductor
# scuttle              * 1.8.4     2023-01-19 [2] Bioconductor
# sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.2.1)
# SingleCellExperiment * 1.20.1    2023-03-17 [2] Bioconductor
# sparseMatrixStats      1.10.0    2022-11-01 [2] Bioconductor
# SummarizedExperiment * 1.28.0    2022-11-01 [2] Bioconductor
# tibble                 3.2.1     2023-03-20 [2] CRAN (R 4.2.3)
# tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.2.1)
# utf8                   1.2.3     2023-01-31 [2] CRAN (R 4.2.2)
# vctrs                  0.6.2     2023-04-19 [1] CRAN (R 4.2.3)
# vipor                  0.4.5     2017-03-22 [2] CRAN (R 4.2.1)
# viridis              * 0.6.2     2021-10-13 [2] CRAN (R 4.2.1) 
