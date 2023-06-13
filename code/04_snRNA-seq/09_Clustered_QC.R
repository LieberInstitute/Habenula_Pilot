## March 2, 2023 - Bukola Ajanaku
# # We've decided (based on heatmaps against Erik's clusters and against our gene
# # marker interest list) to proceed with Walktrap 10 which has 37  clusters.
# # This script will updated the annotations for each cluster based on further investigation
# # and also plot sn-RNA-seq qc against for each cluster to make sure no cluster is
# # driven by poor QC.
# qrsh -l mem_free=20G,h_vmem=20G

library("SingleCellExperiment")
library("here")
library("sessioninfo")
library("xlsx")
library("tidyverse")
library("scater")
library("jaffelab")

# most updaated sce object with first round of annotations (where we selected to move forward with wT10)
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
          "sce_post_clustering_with_logcounts.Rdata"))

# pulling excel sheet with updated annotations (looking for snType and bulkType)
snAnno <- read.xlsx(file = here("processed-data", "04_snRNA-seq", "09_Clustered_QC",
                                  "GeneAnnotations_Bukola.xlsx"), sheetName = c("WalkTrap10"))

# cleaning up for standardization
snAnno_clean <- snAnno |> 
  filter(!is.na(Cluster)) |>
  mutate(Type_clean = gsub("[^a-zA-Z0-9]", "", Type), 
         Cluster = paste0("10wTrap_", Cluster)) |>
  mutate(splitProbClusts = paste(problemClusts, ss(Cluster, "_", 2), sep = "_")) |>
  mutate(splitSNType = paste(snType, ss(Cluster, "_", 2), sep = "_"))

# dropping three extra columns 
snAnno_clean <- select(snAnno_clean, -c(starts_with("Na"), "MoreInfo"))

# sanity check (checking for repeats in annotations to check that split by cluster makes sense)
snAnno_clean |> count(snType)
snAnno_clean |> count(problemClusts)

# using match to add annotated name columns for wt10
## actual annotatios for now  
  sce$snAnno <- snAnno_clean$snType[match(sce$wT_10_Erik, snAnno_clean$Cluster)]

## adding columns that we want to plot (annotations but separated by clusters and 
## annotations with problematic cluster highlights all split by clusters)
  sce$splitProbClusts <- snAnno_clean$splitProbClusts[match(sce$wT_10_Erik, snAnno_clean$Cluster)]
  sce$splitSNType <- snAnno_clean$splitSNType[match(sce$wT_10_Erik, snAnno_clean$Cluster)]
  
# sourcing for custom  color palette 
source(here("code", "04_snRNA-seq", "sourcing", "color_Scheme_CT.R"))

# Now that sce has been established with new annotations from Excel sheet, we can
# plot against qc metrics (also in colData of sce with annotation identities):
####### Exploring Barcodes  #####################################################
# always create plot dir before plotting 
  plot_dir <- here("plots", "04_snRNA-seq", "09_Clustered_QC")
  if(!dir.exists(plot_dir)){
    dir.create(plot_dir)
  } 


# Checking library size ("sum": the higher the better, dropping lows)
pdf(file = here(plot_dir, "wt10_LIBSIZE_QC.pdf"), height = 7, width = 11)
  # coloring by annotations with problems highlighted
  ggcells(sce, mapping = aes(x = splitProbClusts, y = sum, fill = splitProbClusts)) +
      geom_boxplot() +
      scale_fill_manual(values = grabColors(length(unique(sce$splitProbClusts)))) +
      theme_linedraw() +
      theme(
        legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
      ggtitle("Library Size with Problematic Cluster Highlights")
  
  # coloring by annotations only (not split)
  ggcells(sce, mapping = aes(x = snAnno, y = sum, fill = snAnno)) +
    geom_boxplot() +
    scale_fill_manual(values = grabColors(length(unique(sce$snAnno)))) +
    theme_linedraw() +
    theme(
      legend.position = "none", axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
    ggtitle("Library Size for all Annotated Clusters")
  
  # regular annotations split by cluster 
  ggcells(sce, mapping = aes(x = splitSNType, y = sum, fill = splitSNType)) +
    geom_boxplot() +
    scale_fill_manual(values = grabColors(length(unique(sce$splitSNType)))) +
    theme_linedraw() +
    theme(
      legend.position = "none", axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle("Library Size for Regular Annotations Split by Cluster")
dev.off()

# Checking detected features ("detected": the higher the better, dropping lows)
pdf(file = here(plot_dir, "wt10_DETECTED_QC.pdf"), height = 7, width = 11)
  # coloring by annotations with problems highlighted
  ggcells(sce, mapping = aes(x = splitProbClusts, y = detected, fill = splitProbClusts)) +
    geom_boxplot() +
    scale_fill_manual(values = grabColors(length(unique(sce$splitProbClusts)))) +
    theme_linedraw() +
    theme(
      legend.position = "none", axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle("Detected Features with Problematic Cluster Highlights")
  
  # coloring by annotations only (not split)
  ggcells(sce, mapping = aes(x = snAnno, y = detected, fill = snAnno)) +
    geom_boxplot() +
    scale_fill_manual(values = grabColors(length(unique(sce$snAnno)))) +
    theme_linedraw() +
    theme(
      legend.position = "none", axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle("Detected Features for all Annotated Clusters")
  
  # regular annotations split by cluster 
  ggcells(sce, mapping = aes(x = splitSNType, y = detected, fill = splitSNType)) +
    geom_boxplot() +
    scale_fill_manual(values = grabColors(length(unique(sce$splitSNType)))) +
    theme_linedraw() +
    theme(
      legend.position = "none", axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle("Detected Features for Regular Annotations Split by Cluster")
  
dev.off()

# Checking mitochondrial rate ("subsets_Mito_percent": the lower the better, dropping highs)
pdf(file = here(plot_dir, "wt10_MITO_PERCENT_QC.pdf"), height = 7, width = 11)
# coloring by annotations with problems highlighted
    ggcells(sce, mapping = aes(x = splitProbClusts, y = subsets_Mito_percent, fill = splitProbClusts)) +
      geom_boxplot() +
      scale_fill_manual(values = grabColors(length(unique(sce$splitProbClusts)))) +
      theme_linedraw() +
      theme(
        legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      ggtitle("Mito Percent with Problematic Cluster Highlights")
    
    # coloring by annotations only (not split)
    ggcells(sce, mapping = aes(x = snAnno, y = subsets_Mito_percent, fill = snAnno)) +
      geom_boxplot() +
      scale_fill_manual(values = grabColors(length(unique(sce$snAnno)))) +
      theme_linedraw() +
      theme(
        legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      ggtitle("Mito Percent for all Annotated Clusters")
    
    # regular annotations split by cluster 
    ggcells(sce, mapping = aes(x = splitSNType, y = subsets_Mito_percent, fill = splitSNType)) +
      geom_boxplot() +
      scale_fill_manual(values = grabColors(length(unique(sce$splitSNType)))) +
      theme_linedraw() +
      theme(
        legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      ggtitle("Mito Percent for Regular Annotations Split by Cluster")
dev.off()

# Checking doublelt scores ("doubletScore": the lower the better, typically dropping anything over 5)
pdf(file = here(plot_dir, "wt10_DOUBLET_SCORE_QC.pdf"), height = 7, width = 11)
# coloring by annotations with problems highlighted
    ggcells(sce, mapping = aes(x = splitProbClusts, y = doubletScore, fill = splitProbClusts)) +
      geom_boxplot() +
      scale_fill_manual(values = grabColors(length(unique(sce$splitProbClusts)))) +
      theme_linedraw() +
      theme(
        legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      ggtitle("Doublet Score with Problematic Cluster Highlights")
    
    # coloring by annotations only (not split)
    ggcells(sce, mapping = aes(x = snAnno, y = doubletScore, fill = snAnno)) +
      geom_boxplot() +
      scale_fill_manual(values = grabColors(length(unique(sce$snAnno)))) +
      theme_linedraw() +
      theme(
        legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      ggtitle("Doublet Score for all Annotated Clusters")
    
    # regular annotations split by cluster 
    ggcells(sce, mapping = aes(x = splitSNType, y = doubletScore, fill = splitSNType)) +
      geom_boxplot() +
      scale_fill_manual(values = grabColors(length(unique(sce$splitSNType)))) +
      theme_linedraw() +
      theme(
        legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      ggtitle("Doublet Score for Regular Annotations Split by Cluster")
dev.off()


# saving sce object
save(sce, file = here("processed-data", "04_snRNA-seq", "sce_objects", 
                                "sce_post_09_clustered_qc.Rdata"))

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
# DelayedArray           0.24.0    2022-11-01 [2] Bioconductor
# DelayedMatrixStats     1.20.0    2022-11-01 [2] Bioconductor
# dplyr                * 1.1.1     2023-03-22 [2] CRAN (R 4.2.3)
# fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.2.2)
# forcats              * 1.0.0     2023-01-29 [2] CRAN (R 4.2.2)
# fs                     1.6.1     2023-02-06 [2] CRAN (R 4.2.2)
# gargle                 1.3.0     2023-01-30 [2] CRAN (R 4.2.2)
# generics               0.1.3     2022-07-05 [2] CRAN (R 4.2.1)
# GenomeInfoDb         * 1.34.9    2023-02-02 [2] Bioconductor
# GenomeInfoDbData       1.2.9     2022-09-29 [2] Bioconductor
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
# IRanges              * 2.32.0    2022-11-01 [2] Bioconductor
# irlba                  2.3.5.1   2022-10-03 [2] CRAN (R 4.2.1)
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
# purrr                * 1.0.1     2023-01-10 [2] CRAN (R 4.2.2)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.2.1)
# rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 4.2.2)
# RColorBrewer           1.1-3     2022-04-03 [2] CRAN (R 4.2.1)
# Rcpp                   1.0.10    2023-01-22 [2] CRAN (R 4.2.2)
# RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.2.3)
# readr                * 2.1.4     2023-02-10 [2] CRAN (R 4.2.2)
# rJava                  1.0-6     2021-12-10 [2] CRAN (R 4.2.1)
# rlang                  1.1.1     2023-04-28 [1] CRAN (R 4.2.3)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.2.1)
# rsvd                   1.0.5     2021-04-16 [2] CRAN (R 4.2.1)
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
# vipor                  0.4.5     2017-03-22 [2] CRAN (R 4.2.1)
# viridis                0.6.2     2021-10-13 [2] CRAN (R 4.2.1)
# viridisLite            0.4.2     2023-05-02 [1] CRAN (R 4.2.3)
# withr                  2.5.0     2022-03-03 [2] CRAN (R 4.2.1)