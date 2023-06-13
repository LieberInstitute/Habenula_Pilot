# January 17, 2023 - Bukola Ajanaku
# Creating pca plots for filtered (removed empty droplets), quality controlled 
# (dropped high mito, low library size, low detected features, and any genes 
# with 0 counts across samples) sce object.
# Based on:
# 1) https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/code/03_build_sce/04_GLM_PCA.R
# 2) https://www.stephaniehicks.com/biocdemo/articles/Demo.html
# qrsh -l mem_free=50G,h_vmem=50G

library("SingleCellExperiment")
library("scran")
library("scater")
library("scry")
library("here")
library("sessioninfo")
library("mbkmeans")
library("HDF5Array")
library("ggplot2")
library("dplyr")
library("scales")
library("viridis")

# loading post qc sce object 
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
           "sce_hb_postQC.Rdata"))
sce <- sce_hb_postQC
rm(sce_hb_postQC)

####### Adding new metric: Run Data (before subsetting) ########################
# Adding Run data 
colData(sce)$Run <- NA
colData(sce)[colData(sce)$Sample == "Br1469", "Run"] <- 1
colData(sce)[colData(sce)$Sample %in% c("Br5558", "Br1204"), "Run"] <- 2
colData(sce)[colData(sce)$Sample %in% c("Br1092", "Br1735", "Br5555", "Br5639"), "Run"] <- 3
sce$Run <- as.factor(sce$Run)

################################################################################

# Deviance featuring selection
sce <- devianceFeatureSelection(sce,
        assay = "counts", fam = "binomial", sorted = F,
        batch = as.factor(sce$Sample))

# Checking outputs for deviance ft selection function
pdf(here("plots", "04_snRNA-seq", "04_GLM_PCA_plots", "binomial_deviance_uncorrected.pdf"))
  plot(sort(rowData(sce)$binomial_deviance, decreasing = T),
       type = "l", xlab = "ranked genes",
       ylab = "binomial deviance", main = "Feature Selection with Deviance"
  )
  abline(v = 2000, lty = 2, col = "red")
dev.off()

# Taking the GLM_PCA approach
sce <- nullResiduals(sce,
                     assay = "counts", fam = "binomial", # default params
                     type = "deviance")

# Selects for the top 2000 most variable genes
hdgs.hb <- rownames(sce)[order(rowData(sce)$binomial_deviance, decreasing = T)][1:2000]

# RUNNING PCA --
sce_uncorrected <- runPCA(sce,
                          exprs_values = "binomial_deviance_residuals",
                          subset_row = hdgs.hb, ncomponents = 100,
                          name = "GLMPCA_approx",
                          BSPARAM = BiocSingular::IrlbaParam()
)

# Plotting using GLM-PCA functions [color by Sample]
pdf(here("plots", "04_snRNA-seq", "04_GLM_PCA_plots", "GLM_uncorrected_plot_by_Sample.pdf"))
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "Sample") 
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "Sample") + 
    facet_wrap(~ sce_uncorrected$Sample)
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "Sample", ncomponents = c(2,3)) # looks batchy
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "Sample", ncomponents = c(2,3)) + 
    facet_wrap(~ sce_uncorrected$Sample)
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "Sample", ncomponents = 5)
dev.off()

# Plotting using GLM-PCA functions [color by Run]
pdf(here("plots", "04_snRNA-seq", "04_GLM_PCA_plots", "GLM_uncorrected_plot_by_Run.pdf"))
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "Run") 
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "Run") + 
    facet_wrap(~ sce_uncorrected$Sample)              
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "Run", ncomponents = c(2,3)) 
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "Run", ncomponents = c(2,3)) + 
    facet_wrap(~ sce_uncorrected$Sample)
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "Run", ncomponents = 5)
dev.off()

# Plotting using GLM-PCA functions [color by Erik Cluster]
pdf(here("plots", "04_snRNA-seq", "04_GLM_PCA_plots", "GLM_uncorrected_plot_by_ct_Erik.pdf"))
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "ct_Erik") + labs(caption = "887 NAs")
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "ct_Erik") + 
    facet_wrap(~ sce_uncorrected$Sample)  + labs(caption = "887 NAs")            
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "ct_Erik", ncomponents = c(2,3)) + 
    labs(caption = "887 NAs")
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "ct_Erik", ncomponents = c(2,3)) + 
    facet_wrap(~ sce_uncorrected$Sample) + labs(caption = "887 NAs")
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "ct_Erik", ncomponents = 5) + 
    labs(caption = "887 NAs")
dev.off()

# Plotting by continuous data
pdf(here("plots", "04_snRNA-seq", "04_GLM_PCA_plots", "GLM_uncorrected_continous_metrics.pdf"))
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "sum") +
    scale_color_viridis(option = "G", begin = 0.19) 
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "sum")  + 
    facet_wrap(~ sce_uncorrected$Sample) +
    scale_color_viridis(option = "G", begin = 0.19) 
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "detected") +
    scale_color_viridis(option = "G", begin = 0.19) 
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "detected")  + 
    facet_wrap(~ sce_uncorrected$Sample) +
    scale_color_viridis(option = "G", begin = 0.19) 
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "subsets_Mito_percent") +
    scale_color_viridis(option = "G", begin = 0.19) 
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "subsets_Mito_percent") + 
    facet_wrap(~ sce_uncorrected$Sample) +
    scale_color_viridis(option = "G", begin = 0.19) 
dev.off()

# RUNNING TSNE and UMAP --
sce_uncorrected <- runTSNE(sce_uncorrected, dimred = "GLMPCA_approx")
sce_uncorrected <- runUMAP(sce_uncorrected, dimred = "GLMPCA_approx")

# Plotting TSNE 
pdf(here("plots", "04_snRNA-seq", "04_GLM_PCA_plots", "TSNE_uncorrected_plot.pdf"))
  plotReducedDim(sce_uncorrected, dimred = "TSNE", colour_by = "Sample") 
  plotReducedDim(sce_uncorrected, dimred = "TSNE", colour_by = "Sample") + facet_wrap(~ sce_uncorrected$Sample)
  plotReducedDim(sce_uncorrected, dimred = "TSNE", colour_by = "Run") 
  plotReducedDim(sce_uncorrected, dimred = "TSNE", colour_by = "Run") + facet_wrap(~ sce_uncorrected$Sample)
  plotReducedDim(sce_uncorrected, dimred = "TSNE", colour_by = "ct_Erik") + 
    labs(caption = "887 NAs")
  plotReducedDim(sce_uncorrected, dimred = "TSNE", colour_by = "ct_Erik") + 
    facet_wrap(~ sce_uncorrected$Sample) + 
    labs(caption = "887 NAs")
dev.off()

# Plotting UMAP
pdf(here("plots", "04_snRNA-seq", "04_GLM_PCA_plots", "UMAP_uncorrected_plot.pdf"))
  plotReducedDim(sce_uncorrected, dimred = "UMAP", colour_by = "Sample")
  plotReducedDim(sce_uncorrected, dimred = "UMAP", colour_by = "Sample") + facet_wrap(~ sce_uncorrected$Sample)
  plotReducedDim(sce_uncorrected, dimred = "UMAP", colour_by = "Run") 
  plotReducedDim(sce_uncorrected, dimred = "UMAP", colour_by = "Run") + facet_wrap(~ sce_uncorrected$Sample)
  plotReducedDim(sce_uncorrected, dimred = "UMAP", colour_by = "ct_Erik") + 
    labs(caption = "887 NAs")
  plotReducedDim(sce_uncorrected, dimred = "UMAP", colour_by = "ct_Erik") + 
    facet_wrap(~ sce_uncorrected$Sample) + 
    labs(caption = "887 NAs")
dev.off()

# Plotting by continuous data, TSNE
pdf(here("plots", "04_snRNA-seq", "04_GLM_PCA_plots", "TSNE_uncorrected_plot_continous_mets.pdf"))
  plotReducedDim(sce_uncorrected, dimred = "TSNE", colour_by = "sum") +
    scale_color_viridis(option = "G", begin = 0.19) 
  plotReducedDim(sce_uncorrected, dimred = "TSNE", colour_by = "sum") + 
    facet_wrap(~ sce_uncorrected$Sample) +
    scale_color_viridis(option = "G", begin = 0.19) 
  plotReducedDim(sce_uncorrected, dimred = "TSNE", colour_by = "detected") +
    scale_color_viridis(option = "G", begin = 0.19) 
  plotReducedDim(sce_uncorrected, dimred = "TSNE", colour_by = "detected") + 
    facet_wrap(~ sce_uncorrected$Sample) +
    scale_color_viridis(option = "G", begin = 0.19) 
  plotReducedDim(sce_uncorrected, dimred = "TSNE", colour_by = "detected") +
    scale_color_viridis(option = "G", begin = 0.19) 
  plotReducedDim(sce_uncorrected, dimred = "TSNE", colour_by = "subsets_Mito_percent") + 
    facet_wrap(~ sce_uncorrected$Sample) +
    scale_color_viridis(option = "G", begin = 0.19) 
dev.off()

# Plotting by continous data, UMAP
pdf(here("plots", "04_snRNA-seq", "04_GLM_PCA_plots", "UMAP_uncorrected_plot_continous_mets.pdf"))
  plotReducedDim(sce_uncorrected, dimred = "UMAP", colour_by = "sum") +
    scale_color_viridis(option = "G", begin = 0.19) 
  plotReducedDim(sce_uncorrected, dimred = "UMAP", colour_by = "sum") + 
    facet_wrap(~ sce_uncorrected$Sample) +
    scale_color_viridis(option = "G", begin = 0.19) 
  plotReducedDim(sce_uncorrected, dimred = "UMAP", colour_by = "detected") +
    scale_color_viridis(option = "G", begin = 0.19)  
  plotReducedDim(sce_uncorrected, dimred = "UMAP", colour_by = "detected") + 
    facet_wrap(~ sce_uncorrected$Sample) +
    scale_color_viridis(option = "G", begin = 0.19) 
  plotReducedDim(sce_uncorrected, dimred = "UMAP", colour_by = "subsets_Mito_percent")+
    scale_color_viridis(option = "G", begin = 0.19)  
  plotReducedDim(sce_uncorrected, dimred = "UMAP", colour_by = "subsets_Mito_percent") +
    facet_wrap(~ sce_uncorrected$Sample) +
    scale_color_viridis(option = "G", begin = 0.19) 
dev.off()

## Saving uncorrected sce object post pca 
save(sce_uncorrected, file = here("processed-data", "04_snRNA-seq", "sce_objects", 
                                  "sce_uncorrected_PCA.Rdata"))

# Saving as HDF5 for later clustering
saveHDF5SummarizedExperiment(sce_uncorrected, dir = here("processed-data", "04_snRNA-seq", 
                            "sce_objects", "sce_uncorrected"), replace = TRUE)

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
# benchmarkme            1.0.8     2022-06-12 [2] CRAN (R 4.2.1)
# benchmarkmeData        1.0.4     2020-04-23 [2] CRAN (R 4.2.1)
# Biobase              * 2.58.0    2022-11-01 [2] Bioconductor
# BiocGenerics         * 0.44.0    2022-11-01 [2] Bioconductor
# BiocNeighbors          1.16.0    2022-11-01 [2] Bioconductor
# BiocParallel           1.32.6    2023-03-17 [2] Bioconductor
# BiocSingular           1.14.0    2022-11-01 [2] Bioconductor
# bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.2.1)
# bluster                1.8.0     2022-11-01 [2] Bioconductor
# cli                    3.6.1     2023-03-23 [2] CRAN (R 4.2.3)
# cluster                2.1.4     2022-08-22 [3] CRAN (R 4.2.3)
# ClusterR               1.3.0     2023-01-21 [2] CRAN (R 4.2.2)
# codetools              0.2-19    2023-02-01 [3] CRAN (R 4.2.3)
# colorout             * 1.2-2     2023-02-09 [1] Github (jalvesaq/colorout@79931fd)
# colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.2.2)
# DelayedArray         * 0.24.0    2022-11-01 [2] Bioconductor
# DelayedMatrixStats     1.20.0    2022-11-01 [2] Bioconductor
# doParallel             1.0.17    2022-02-07 [2] CRAN (R 4.2.1)
# dplyr                * 1.1.1     2023-03-22 [2] CRAN (R 4.2.3)
# dqrng                  0.3.0     2021-05-01 [2] CRAN (R 4.2.1)
# edgeR                  3.40.2    2023-01-19 [2] Bioconductor
# fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.2.2)
# foreach                1.5.2     2022-02-02 [2] CRAN (R 4.2.1)
# generics               0.1.3     2022-07-05 [2] CRAN (R 4.2.1)
# GenomeInfoDb         * 1.34.9    2023-02-02 [2] Bioconductor
# GenomeInfoDbData       1.2.9     2022-09-29 [2] Bioconductor
# GenomicRanges        * 1.50.2    2022-12-16 [2] Bioconductor
# ggbeeswarm             0.7.1     2022-12-16 [2] CRAN (R 4.2.2)
# ggplot2              * 3.4.2     2023-04-03 [2] CRAN (R 4.2.3)
# ggrepel                0.9.3     2023-02-03 [2] CRAN (R 4.2.2)
# glue                   1.6.2     2022-02-24 [2] CRAN (R 4.2.1)
# gmp                    0.7-1     2023-02-07 [2] CRAN (R 4.2.2)
# gridExtra              2.3       2017-09-09 [2] CRAN (R 4.2.1)
# gtable                 0.3.3     2023-03-21 [2] CRAN (R 4.2.3)
# HDF5Array            * 1.26.0    2022-11-01 [2] Bioconductor
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.2.1)
# httr                   1.4.5     2023-02-24 [2] CRAN (R 4.2.2)
# igraph                 1.4.2     2023-04-07 [2] CRAN (R 4.2.3)
# IRanges              * 2.32.0    2022-11-01 [2] Bioconductor
# irlba                  2.3.5.1   2022-10-03 [2] CRAN (R 4.2.1)
# iterators              1.0.14    2022-02-05 [2] CRAN (R 4.2.1)
# lattice                0.20-45   2021-09-22 [3] CRAN (R 4.2.3)
# lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.2.1)
# limma                  3.54.2    2023-02-28 [2] Bioconductor
# locfit                 1.5-9.7   2023-01-02 [2] CRAN (R 4.2.2)
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.2.1)
# Matrix               * 1.5-4     2023-04-04 [2] CRAN (R 4.2.3)
# MatrixGenerics       * 1.10.0    2022-11-01 [2] Bioconductor
# matrixStats          * 1.0.0     2023-06-02 [1] CRAN (R 4.2.3)
# mbkmeans             * 1.14.0    2022-11-01 [2] Bioconductor
# metapod                1.6.0     2022-11-01 [2] Bioconductor
# munsell                0.5.0     2018-06-12 [2] CRAN (R 4.2.1)
# pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.2.3)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.2.1)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.2.1)
# Rcpp                   1.0.10    2023-01-22 [2] CRAN (R 4.2.2)
# RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.2.3)
# rhdf5                * 2.42.1    2023-04-07 [2] Bioconductor
# rhdf5filters           1.10.1    2023-03-24 [2] Bioconductor
# Rhdf5lib               1.20.0    2022-11-01 [2] Bioconductor
# rlang                  1.1.1     2023-04-28 [1] CRAN (R 4.2.3)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.2.1)
# rsvd                   1.0.5     2021-04-16 [2] CRAN (R 4.2.1)
# S4Vectors            * 0.36.2    2023-02-26 [2] Bioconductor
# ScaledMatrix           1.6.0     2022-11-01 [2] Bioconductor
# scales               * 1.2.1     2022-08-20 [2] CRAN (R 4.2.1)
# scater               * 1.26.1    2022-11-13 [2] Bioconductor
# scran                * 1.26.2    2023-01-19 [2] Bioconductor
# scry                 * 1.10.0    2022-11-01 [2] Bioconductor
# scuttle              * 1.8.4     2023-01-19 [2] Bioconductor
# sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.2.1)
# SingleCellExperiment * 1.20.1    2023-03-17 [2] Bioconductor
# sparseMatrixStats      1.10.0    2022-11-01 [2] Bioconductor
# statmod                1.5.0     2023-01-06 [2] CRAN (R 4.2.2)
# SummarizedExperiment * 1.28.0    2022-11-01 [2] Bioconductor
# tibble                 3.2.1     2023-03-20 [2] CRAN (R 4.2.3)
# tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.2.1)
# utf8                   1.2.3     2023-01-31 [2] CRAN (R 4.2.2)
# vctrs                  0.6.2     2023-04-19 [1] CRAN (R 4.2.3)
# vipor                  0.4.5     2017-03-22 [2] CRAN (R 4.2.1)
# viridis              * 0.6.2     2021-10-13 [2] CRAN (R 4.2.1)
# viridisLite          * 0.4.2     2023-05-02 [1] CRAN (R 4.2.3)
# withr                  2.5.0     2022-03-03 [2] CRAN (R 4.2.1)
# XVector                0.38.0    2022-11-01 [2] Bioconductor
# zlibbioc               1.44.0    2022-11-01 [2] Bioconductor
# 
# [1] /users/bsimbiat/R/4.2.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/library
