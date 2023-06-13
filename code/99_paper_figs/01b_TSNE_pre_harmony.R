## April 17, 2023 - Bukola Ajanaku
# Plotting same TSNE plots but pre-harmony.
# qrsh -l mem_free=50G,h_vmem=50G

library("SingleCellExperiment")
library("here")
library("tidyverse")
library("scater")
library("ggplot2")
library("cowplot")

# loading old sce object (post qc sce object)
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
         "sce_uncorrected_PCA.Rdata"), verbose = TRUE)
# sce_uncorrected 

dim(sce_uncorrected)
# [1] 33848 17082

# loading official sce object 
load(here("processed-data", "99_paper_figs", "sce_objects", 
          "official_final_sce.RDATA"), verbose = TRUE)
# sce

dim(sce)
# [1] 33848 17031

# making sure colnames of sce_uncorrected are unique 
colnames(sce_uncorrected) <- paste0(sce_uncorrected$Sample, "_", sce_uncorrected$Barcode)

# subsetting sce_uncorrected to only the nuclei we've kept in sce
sce_uncorrected_clean <- sce_uncorrected[, which(colnames(sce_uncorrected) %in% colnames(sce))]

dim(sce_uncorrected_clean)
# [1] 33848 17031

# Now making sure phenotype data from sce is in sce_uncorrected (clean)
colData(sce_uncorrected_clean) <- colData(sce)

# just checking for final Annotations
table(sce_uncorrected_clean$final_Annotations)
    # Astrocyte       Endo Excit.Thal Inhib.Thal      LHb.1      LHb.2      LHb.3 
    # 538         38       1800       7612        201        266        134 
    # LHb.4      LHb.5      LHb.6      LHb.7      MHb.1      MHb.2      MHb.3 
    # 477         83         39       1014        152        540         18 
    # Microglia      Oligo        OPC 
    # 145       2178       1796 

#### Prepping for plotting #####################################################
# creating plot directory
plot_dir <- here("plots", "99_paper_figs", "01_TSNEs", "Pre-Harmony")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

# adding color pallete (same color scheme used for progress report heatmap)
cluster_colors <- c( "Oligo" = c("#475c6c"), 
                     "OPC"= c("#899499"), 
                     "Microglia" = c("#bfb5b2"), 
                     "Astrocyte" = c("#c7bbc9"), 
                     "Endo" = c("#8a8583"), 
                     "Excit.Neuron" = c("#cd8b62"), 
                     "Inhib.Thal" = c("#eed7a1"),  
                     "Excit.Thal" = c('#f7efd2'), 
                     "LHb.1" = c("#00FFFF"),
                     "LHb.2" = c("#0096FF"), 
                     "LHb.3" = c ("#1434A4"), 
                     "LHb.4" = c("#00008B"), 
                     "LHb.5" = c("#40E0D0"), 
                     "LHb.6" = c("#008080"),  
                     "LHb.7" = c("#7DF9FF"), 
                     "MHb.1" = c("#800020"), 
                     "MHb.2" = c("#D70040"),
                     "MHb.3" = c("#D2042D") 
)

## create sce_scorted and unsorted based on NeuN
sce_unc_sorted <- sce_uncorrected_clean[, which(sce_uncorrected_clean$NeuN == "NeuN.Sorted")]
sce_unc_unsorted <- sce_uncorrected_clean[, which(sce_uncorrected_clean$NeuN == "NeuN.Unsorted")]

##### PLOTTING TSNEs ###########################################################
# Pre-Harmonization colored by Samplel and faceted by Sample
pdf(here(plot_dir, "TSNE_uncorrected_by_Sample.pdf"), width = 14, height = 9)
plot1 <-  plotReducedDim(sce_uncorrected_clean, dimred = "TSNE") +
    geom_point(aes(color = sce_uncorrected_clean$Sample)) + 
    guides(color = guide_legend(title="Sample ID")) + 
  theme(legend.position = "none")

plot2 <-  plotReducedDim(sce_uncorrected_clean, dimred = "TSNE") +
  geom_point(aes(color = sce_uncorrected_clean$Sample)) + 
  guides(color = guide_legend(title="Sample ID")) +
  facet_wrap(~ sce_uncorrected_clean$Sample)

plot_grid(plot1, plot2)
dev.off()

# Pre-Harmonization colored by Samplel and faceted by Run
pdf(here(plot_dir, "TSNE_uncorrected_by_Run.pdf"), 
    width = 14, height = 6)
plot1 <- plotReducedDim(sce_uncorrected_clean, dimred = "TSNE") +
  geom_point(aes(color = sce_uncorrected_clean$Sample)) + 
  theme(legend.position = "none")

plot2 <- plotReducedDim(sce_uncorrected_clean, dimred = "TSNE") +
  geom_point(aes(color = sce_uncorrected_clean$Sample)) + 
  facet_wrap(~ sce_uncorrected_clean$Run) + 
  guides(color = guide_legend(title="Sample ID"))

plot_grid(plot1, plot2, rel_widths = c(1,2))
dev.off()

# for NeuN sorting and Non-NeuN sorting
# plot_sorted <- plotReducedDim(sce_unc_sorted, dimred = "TSNE") +
#   geom_point(aes(color = sce_unc_sorted$final_Annotations), alpha = 0.4) + 
#   scale_colour_manual(values = cluster_colors) +
#   facet_grid(sce_unc_sorted$NeuN ~ sce_unc_sorted$Sample) + 
#   guides(color = guide_legend(title="Cell Type"))
# 
# plot_unsorted <-   plotReducedDim(sce_unc_unsorted, dimred = "TSNE") +
#   geom_point(aes(color = sce_unc_unsorted$final_Annotations), alpha = 0.4) + 
#   scale_colour_manual(values = cluster_colors) +
#   facet_grid(sce_unc_unsorted$NeuN ~ sce_unc_unsorted$Sample) + 
#   guides(color = guide_legend(title="Cell Type"))
# 
# 
# pdf(here(plot_dir, "TSNE_harmony_by_finalAnno_splitbySampleAndSorting_PRE-HARMONY.pdf"), width = 13, height = 9)
# plot_grid(
#   plot_sorted,
#   plot_unsorted,
#   ncol = 1
# )
# dev.off()
# 
# pdf(here(plot_dir, "TSNE_harmony_by_final_Annotations_splitbyRun_PRE-HARMONY.pdf"), 
#     width = 20, height = 8)
# plot1 <- plotReducedDim(sce_uncorrected_clean, dimred = "TSNE") +
#   geom_point(aes(color = sce_uncorrected_clean$final_Annotations)) + 
#   scale_colour_manual(values = cluster_colors) +
#   facet_wrap(~ sce_uncorrected_clean$Run) + 
#   theme(legend.position = "none")
# 
# plot2 <- plotReducedDim(sce_uncorrected_clean, dimred = "TSNE") +
#   geom_point(aes(color = sce_uncorrected_clean$final_Annotations)) + 
#   scale_colour_manual(values = cluster_colors) +
#   facet_wrap(~ sce_uncorrected_clean$final_Annotations) + 
#   guides(color = guide_legend(title="Cell Type"))
# 
# plot_grid(plot1, plot2)
# dev.off()

# Saving uncorrected sce object 
save(sce_uncorrected_clean, file = here("processed-data", "99_paper_figs", "sce_objects", 
           "official_final_uncorrected_sce.RDATA"))

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
# dplyr                * 1.1.1     2023-03-22 [2] CRAN (R 4.2.3)
# fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.2.2)
# forcats              * 1.0.0     2023-01-29 [2] CRAN (R 4.2.2)
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
# hms                    1.1.3     2023-03-21 [2] CRAN (R 4.2.3)
# IRanges              * 2.32.0    2022-11-01 [2] Bioconductor
# irlba                  2.3.5.1   2022-10-03 [2] CRAN (R 4.2.1)
# lattice                0.20-45   2021-09-22 [3] CRAN (R 4.2.3)
# lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.2.1)
# lubridate            * 1.9.2     2023-02-10 [2] CRAN (R 4.2.2)
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.2.1)
# Matrix                 1.5-4     2023-04-04 [2] CRAN (R 4.2.3)
# MatrixGenerics       * 1.10.0    2022-11-01 [2] Bioconductor
# matrixStats          * 1.0.0     2023-06-02 [1] CRAN (R 4.2.3)
# munsell                0.5.0     2018-06-12 [2] CRAN (R 4.2.1)
# pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.2.3)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.2.1)
# purrr                * 1.0.1     2023-01-10 [2] CRAN (R 4.2.2)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.2.1)
# Rcpp                   1.0.10    2023-01-22 [2] CRAN (R 4.2.2)
# RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.2.3)
# readr                * 2.1.4     2023-02-10 [2] CRAN (R 4.2.2)
# rlang                  1.1.1     2023-04-28 [1] CRAN (R 4.2.3)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.2.1)
# rsvd                   1.0.5     2021-04-16 [2] CRAN (R 4.2.1)
# S4Vectors            * 0.36.2    2023-02-26 [2] Bioconductor
# ScaledMatrix           1.6.0     2022-11-01 [2] Bioconductor
# scales                 1.2.1     2022-08-20 [2] CRAN (R 4.2.1)
# scater               * 1.26.1    2022-11-13 [2] Bioconductor
# scuttle              * 1.8.4     2023-01-19 [2] Bioconductor
# sessioninfo            1.2.2     2021-12-06 [2] CRAN (R 4.2.1)
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
# XVector                0.38.0    2022-11-01 [2] Bioconductor
# zlibbioc               1.44.0    2022-11-01 [2] Bioconductor
# 
# [1] /users/bsimbiat/R/4.2.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/library
