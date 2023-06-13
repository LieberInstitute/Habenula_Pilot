## May 2, 2023 - Bukola Ajanaku
# Post QC sce object being shown per annotated cluster.
# qrsh -l mem_free=30G,h_vmem=30G

library("SingleCellExperiment")
library("here")
library("ggplot2")
library("sessioninfo")
library("scater")
library("cowplot")

# loading sce object post-cleaning and annotations (includes Hb cluster)
load(here("processed-data", "99_paper_figs",  "sce_objects", "sce_final_preHbdrop.RDATA"),
     verbose = TRUE)
# sce 

# sourcing official color palette 
source(file = here("code", "99_paper_figs", "source_colors.R"))
# bulk_colors and sn_colors

# creating plot_dir
plot_dir <- here("plots", "99_paper_figs", "05_QC_Boxplots")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

# plotting per metric on final_Annotations 
# also has OPC_noisy class
pdf(file = here(plot_dir, "sce_Annotated_with_Ambig_QC_boxplot.pdf"), width = 12, height = 12)

plot1 <- ggcells(sce, mapping = aes(x = final_Annotations, y = doubletScore, fill = final_Annotations)) +
      geom_boxplot() +
      scale_fill_manual(values = sn_colors) +
      theme_linedraw() +
      theme(
        legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) + labs (y = "Doublet Score")

plot2 <- ggcells(sce, mapping = aes(x = final_Annotations, y = subsets_Mito_percent, fill = final_Annotations)) +
  geom_boxplot() +
  scale_fill_manual(values = sn_colors) +
  theme_linedraw() +
  theme(
    legend.position = "none", axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) + labs (y = "Percent Mito")

plot3 <- ggcells(sce, mapping = aes(x = final_Annotations, y = sum, fill = final_Annotations)) +
  geom_boxplot() +
  scale_fill_manual(values = sn_colors) +
  theme_linedraw() +
  theme(
    legend.position = "none", axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) + labs (y = "Library Size")

plot4 <- ggcells(sce, mapping = aes(x = final_Annotations, y = detected, fill = final_Annotations)) +
  geom_boxplot() +
  scale_fill_manual(values = sn_colors) +
  theme_linedraw() +
  theme(
    legend.position = "none", axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) + labs (y = "Detected Features")

plot_grid(plot1, plot2, plot3, plot4,
          labels = c("A", "B", "C", "D"),
          ncol = 1)

dev.off()

# dropping Excit.Neuron
sce_drop <- sce[, which(!sce$final_Annotations == "Excit.Neuron")]
sce_drop <- sce_drop[, which(!sce_drop$final_Annotations == "OPC_noisy")]

# check 
table(sce_drop$final_Annotations)


# plotting per metric on final_Annotations without ambig class
# also has OPC_noisy class
pdf(file = here(plot_dir, "sce_No_Ambig_Annotated_QC_boxplot.pdf"), width = 12, height = 12)

plot1 <- ggcells(sce_drop, mapping = aes(x = final_Annotations, y = doubletScore, fill = final_Annotations)) +
  geom_boxplot() +
  scale_fill_manual(values = sn_colors) +
  theme_linedraw() +
  theme(
    legend.position = "none", axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) + labs (y = "Doublet Score")

plot2 <- ggcells(sce_drop, mapping = aes(x = final_Annotations, y = subsets_Mito_percent, fill = final_Annotations)) +
  geom_boxplot() +
  scale_fill_manual(values = sn_colors) +
  theme_linedraw() +
  theme(
    legend.position = "none", axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) + labs (y = "Percent Mito")

plot3 <- ggcells(sce_drop, mapping = aes(x = final_Annotations, y = sum, fill = final_Annotations)) +
  geom_boxplot() +
  scale_fill_manual(values = sn_colors) +
  theme_linedraw() +
  theme(
    legend.position = "none", axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) + labs (y = "Library Size")

plot4 <- ggcells(sce_drop, mapping = aes(x = final_Annotations, y = detected, fill = final_Annotations)) +
  geom_boxplot() +
  scale_fill_manual(values = sn_colors) +
  theme_linedraw() +
  theme(
    legend.position = "none", axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) + labs (y = "Detected Features")

plot_grid(plot1, plot2, plot3, plot4,
          labels = c("A", "B", "C", "D"),
          ncol = 1)

dev.off()

##### FOR ONE DRIVE PIECES #####################################################
# dirty qc 
pdf(file = here(plot_dir, "forOneDrive", "sfigu_doublet_score_boxplot_PRE-QC.pdf"), 
    width = 8, height = 3.5)

    ggcells(sce, mapping = aes(x = final_Annotations, y = doubletScore, fill = final_Annotations)) +
      geom_boxplot() +
      scale_fill_manual(values = sn_colors) +
      theme_linedraw() +
      theme(
        legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) + 
      labs (y = "Doublet Score")

dev.off()

pdf(file = here(plot_dir, "forOneDrive", "sfigu_mito_percent_boxplot_PRE-QC.pdf"), 
    width = 8, height = 3.5)

    ggcells(sce, mapping = aes(x = final_Annotations, y = subsets_Mito_percent, fill = final_Annotations)) +
      geom_boxplot() +
      scale_fill_manual(values = sn_colors) +
      theme_linedraw() +
      theme(
        legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) + 
      labs (y = "Percent Mito")

dev.off()

pdf(file = here(plot_dir, "forOneDrive", "sfigu_library_size_boxplot_PRE-QC.pdf"), 
    width = 8, height = 3.5)

    ggcells(sce, mapping = aes(x = final_Annotations, y = sum, fill = final_Annotations)) +
      geom_boxplot() +
      scale_fill_manual(values = sn_colors) +
      theme_linedraw() +
      theme(
        legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) + 
      labs (y = "Library Size")

dev.off()

pdf(file = here(plot_dir, "forOneDrive", "sfigu_detected_features_boxplot_PRE-QC.pdf"), 
    width = 8, height = 3.5)

    ggcells(sce, mapping = aes(x = final_Annotations, y = detected, fill = final_Annotations)) +
      geom_boxplot() +
      scale_fill_manual(values = sn_colors) +
      theme_linedraw() +
      theme(
        legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) + 
      labs (y = "Detected Features")

dev.off()

# post-qc
pdf(file = here(plot_dir, "forOneDrive", "sfigu_doublet_score_boxplot_POST-QCDrop.pdf"), 
    width = 8, height = 3.5)

    ggcells(sce_drop, mapping = aes(x = final_Annotations, y = doubletScore, fill = final_Annotations)) +
      geom_boxplot() +
      scale_fill_manual(values = sn_colors) +
      theme_linedraw() +
      theme(
        legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) + 
      labs (y = "Doublet Score")

dev.off()

pdf(file = here(plot_dir, "forOneDrive", "sfigu_mito_percent_boxplot_POST-QCDrop.pdf"), 
    width = 8, height = 3.5)

    ggcells(sce_drop, mapping = aes(x = final_Annotations, y = subsets_Mito_percent, fill = final_Annotations)) +
      geom_boxplot() +
      scale_fill_manual(values = sn_colors) +
      theme_linedraw() +
      theme(
        legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) + 
      labs (y = "Percent Mito")

dev.off()

pdf(file = here(plot_dir, "forOneDrive", "sfigu_library_size_boxplot_POST-QCDrop.pdf"), 
    width = 8, height = 3.5)

    ggcells(sce_drop, mapping = aes(x = final_Annotations, y = sum, fill = final_Annotations)) +
      geom_boxplot() +
      scale_fill_manual(values = sn_colors) +
      theme_linedraw() +
      theme(
        legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) + 
      labs (y = "Library Size")

dev.off()

pdf(file = here(plot_dir, "forOneDrive", "sfigu_dectected_features_boxplot_POST-QCDrop.pdf"), 
    width = 8, height = 3.5)

    ggcells(sce_drop, mapping = aes(x = final_Annotations, y = detected, fill = final_Annotations)) +
      geom_boxplot() +
      scale_fill_manual(values = sn_colors) +
      theme_linedraw() +
      theme(
        legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) + 
      labs (y = "Detected Features")

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
# cowplot              * 1.1.1     2020-12-30 [2] CRAN (R 4.2.1)
# DelayedArray           0.24.0    2022-11-01 [2] Bioconductor
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
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.2.1)
# IRanges              * 2.32.0    2022-11-01 [2] Bioconductor
# irlba                  2.3.5.1   2022-10-03 [2] CRAN (R 4.2.1)
# lattice                0.20-45   2021-09-22 [3] CRAN (R 4.2.3)
# lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.2.1)
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.2.1)
# Matrix                 1.5-4     2023-04-04 [2] CRAN (R 4.2.3)
# MatrixGenerics       * 1.10.0    2022-11-01 [2] Bioconductor
# matrixStats          * 1.0.0     2023-06-02 [1] CRAN (R 4.2.3)
# munsell                0.5.0     2018-06-12 [2] CRAN (R 4.2.1)
# pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.2.3)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.2.1)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.2.1)
# Rcpp                   1.0.10    2023-01-22 [2] CRAN (R 4.2.2)
# RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.2.3)
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
# viridis                0.6.2     2021-10-13 [2] CRAN (R 4.2.1)
# viridisLite            0.4.2     2023-05-02 [1] CRAN (R 4.2.3)
# withr                  2.5.0     2022-03-03 [2] CRAN (R 4.2.1)
# XVector                0.38.0    2022-11-01 [2] Bioconductor
# zlibbioc               1.44.0    2022-11-01 [2] Bioconductor
# 
# [1] /users/bsimbiat/R/4.2.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/library
# 
# ─────────────────────────────────────────────────────────────────────────────────────────────
