## March 10, 2023 - Bukola Ajanaku
# On snAnno2, I will be running mean expression data before and after combining MHb.2 
# and MHb.3. I will not be renaming any clusters other than those that are being changed in order
# to keep organization clear.
# qrsh -l mem_free=20G,h_vmem=20G

library("SingleCellExperiment")
library("here")
library("tidyverse")
library("sessioninfo")
library("DeconvoBuddies")

# supposed to load from here, but hasn't been saved yet
# load(here("processed-data", "04_snRNA-seq", "sce_objects", "sce_with_snAnno2.RDATA")))

# loading and fixing 
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
          "sce_post_09_clustered_qc.Rdata"))

# creating snAnno2 which will contain the combined MHb groups
sce$snAnno2 <- sce$snAnno

## LHb.6 is actually Endothelial. Total LHb is now 7 from 8.
sce$snAnno2[sce$snAnno2 == "LHb.6"] <- "Endo"

## creating plot_dir for all downstream plots
# for saving plots
plot_dir <- here("plots", "05_explore_sce", "03_Plot_Investigate_snAnno",
                 "05_Updated_Annotations_meanExpression")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

#### BEFORE COMBINING MHB.2 AND MHB.3 ##########################################
# Terminal 4
table(sce$snAnno2)
  # Astrocyte       Endo Excit.Thal         Hb Inhib.Thal      LHb.1      LHb.2 
  # 538         38       1800         51       7612        201        266 
  # LHb.3      LHb.4      LHb.5      LHb.7      LHb.8      MHb.1      MHb.2 
  # 134        477         83         39       1014        152        145 
  # MHb.3      MHb.4  Microglia      Oligo        OPC 
  # 395         18        145       2178       1796 

# grabbing mean ration data
mean_ratio_snAnno2 <- get_mean_ratio2(
  sce,
  cellType_col = "snAnno2",
  assay_name = "logcounts",
  add_symbol = TRUE
)

findMark_manAnnno2 <- findMarkers_1vAll(
  sce,
  assay_name = "counts",
  cellType_col = "snAnno2",
  add_symbol = FALSE,
  mod = "~Sample"
)

## For combined manually annotated clusters 
snAnno_marker_stats_new_2 <- left_join(mean_ratio_snAnno2, findMark_manAnnno2, by = c("gene", "cellType.target"))

pdf(file = here(plot_dir, "Mean_Ratio_Expression_for_new_snAnno2.pdf"))
for (j in levels(as.factor(sce$snAnno2))) {
  message(j)
  
  print(plot_marker_express(sce, 
                             snAnno_marker_stats_new_2, 
                             n_genes = 5,
                             rank_col = "rank_ratio", 
                             anno_col = "anno_ratio",
                             cellType_col = "snAnno2",
                             cell_type = j)
  )
}
dev.off()


#### AFTER COMBINING MHB.2 AND MHB.3 ###########################################
sce$snAnno3 <- sce$snAnno2

# combining MHb.3 with MHb.2
sce$snAnno3[sce$snAnno3 == "MHb.3"] <- "MHb.2"

# RENAMING!
sce$snAnno3[sce$snAnno3 == "LHb.7"] <- "LHb.6"
sce$snAnno3[sce$snAnno3 == "LHb.8"] <- "LHb.7"
sce$snAnno3[sce$snAnno3 == "MHb.4"] <- "MHb.3"

# dropping Hb cluster
sce <- sce[ , which(sce$snAnno3 != "Hb")]

# sanity check
table(sce$snAnno3)

# grabbing mean ration data
mean_ratio_snAnno3_combined <- get_mean_ratio2(
  sce,
  cellType_col = "snAnno3",
  assay_name = "logcounts",
  add_symbol = TRUE
)

findMark_manAnnno3_combined <- findMarkers_1vAll(
  sce,
  assay_name = "counts",
  cellType_col = "snAnno3",
  add_symbol = FALSE,
  mod = "~Sample"
)

## For combined manually annotated clusters 
snAnno_marker_stats_new_3_combined <- left_join(mean_ratio_snAnno3_combined, findMark_manAnnno3_combined, by = c("gene", "cellType.target"))

pdf(file = here(plot_dir, "Mean_Ratio_Expression_for_new_snAnno3_combined.pdf"))
for (j in levels(as.factor(sce$snAnno3))) {
  message(j)
  
  print(plot_marker_express(sce, 
                            snAnno_marker_stats_new_3_combined, 
                            n_genes = 5,
                            rank_col = "rank_ratio", 
                            anno_col = "anno_ratio",
                            cellType_col = "snAnno3",
                            cell_type = j)
  )
}
dev.off()

# saving
save(sce, here("processed-data", "04_snRNA-seq", "sce_objects", "sce_with_snAnno2_and_snAnno3.RDATA"))

save(mean_ratio_snAnno2, findMark_manAnnno2, snAnno_marker_stats_new_2, file =
       here("processed-data", "05_explore_sce", "mean_ratio_for_snAnno2_from_05_Updated_Annotations_meanExpression.Rdata"))
     
save(mean_ratio_snAnno3_combined, findMark_manAnnno3_combined, snAnno_marker_stats_new_3_combined, file =
       here("processed-data", "05_explore_sce", "mean_ratio_for_snAnno3_from_05_Updated_Annotations_meanExpression.Rdata"))

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
# Biobase              * 2.58.0    2022-11-01 [2] Bioconductor
# BiocGenerics         * 0.44.0    2022-11-01 [2] Bioconductor
# BiocNeighbors          1.16.0    2022-11-01 [2] Bioconductor
# BiocParallel           1.32.6    2023-03-17 [2] Bioconductor
# BiocSingular           1.14.0    2022-11-01 [2] Bioconductor
# bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.2.1)
# bluster                1.8.0     2022-11-01 [2] Bioconductor
# cli                    3.6.1     2023-03-23 [2] CRAN (R 4.2.3)
# cluster                2.1.4     2022-08-22 [3] CRAN (R 4.2.3)
# codetools              0.2-19    2023-02-01 [3] CRAN (R 4.2.3)
# colorout             * 1.2-2     2023-02-09 [1] Github (jalvesaq/colorout@79931fd)
# colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.2.2)
# DeconvoBuddies       * 0.99.0    2023-03-10 [1] Github (lahuuki/DeconvoBuddies@d5774e3)
# DelayedArray           0.24.0    2022-11-01 [2] Bioconductor
# DelayedMatrixStats     1.20.0    2022-11-01 [2] Bioconductor
# dplyr                * 1.1.1     2023-03-22 [2] CRAN (R 4.2.3)
# dqrng                  0.3.0     2021-05-01 [2] CRAN (R 4.2.1)
# edgeR                  3.40.2    2023-01-19 [2] Bioconductor
# fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.2.2)
# forcats              * 1.0.0     2023-01-29 [2] CRAN (R 4.2.2)
# generics               0.1.3     2022-07-05 [2] CRAN (R 4.2.1)
# GenomeInfoDb         * 1.34.9    2023-02-02 [2] Bioconductor
# GenomeInfoDbData       1.2.9     2022-09-29 [2] Bioconductor
# GenomicRanges        * 1.50.2    2022-12-16 [2] Bioconductor
# ggplot2              * 3.4.2     2023-04-03 [2] CRAN (R 4.2.3)
# glue                   1.6.2     2022-02-24 [2] CRAN (R 4.2.1)
# gtable                 0.3.3     2023-03-21 [2] CRAN (R 4.2.3)
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.2.1)
# hms                    1.1.3     2023-03-21 [2] CRAN (R 4.2.3)
# igraph                 1.4.2     2023-04-07 [2] CRAN (R 4.2.3)
# IRanges              * 2.32.0    2022-11-01 [2] Bioconductor
# irlba                  2.3.5.1   2022-10-03 [2] CRAN (R 4.2.1)
# lattice                0.20-45   2021-09-22 [3] CRAN (R 4.2.3)
# lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.2.1)
# limma                  3.54.2    2023-02-28 [2] Bioconductor
# locfit                 1.5-9.7   2023-01-02 [2] CRAN (R 4.2.2)
# lubridate            * 1.9.2     2023-02-10 [2] CRAN (R 4.2.2)
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.2.1)
# Matrix                 1.5-4     2023-04-04 [2] CRAN (R 4.2.3)
# MatrixGenerics       * 1.10.0    2022-11-01 [2] Bioconductor
# matrixStats          * 1.0.0     2023-06-02 [1] CRAN (R 4.2.3)
# metapod                1.6.0     2022-11-01 [2] Bioconductor
# munsell                0.5.0     2018-06-12 [2] CRAN (R 4.2.1)
# pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.2.3)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.2.1)
# purrr                * 1.0.1     2023-01-10 [2] CRAN (R 4.2.2)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.2.1)
# rafalib                1.0.0     2015-08-09 [1] CRAN (R 4.2.2)
# RColorBrewer           1.1-3     2022-04-03 [2] CRAN (R 4.2.1)
# Rcpp                   1.0.10    2023-01-22 [2] CRAN (R 4.2.2)
# RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.2.3)
# readr                * 2.1.4     2023-02-10 [2] CRAN (R 4.2.2)
# rlang                  1.1.1     2023-04-28 [1] CRAN (R 4.2.3)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.2.1)
# rsvd                   1.0.5     2021-04-16 [2] CRAN (R 4.2.1)
# S4Vectors            * 0.36.2    2023-02-26 [2] Bioconductor
# ScaledMatrix           1.6.0     2022-11-01 [2] Bioconductor
# scales                 1.2.1     2022-08-20 [2] CRAN (R 4.2.1)
# scran                  1.26.2    2023-01-19 [2] Bioconductor
# scuttle                1.8.4     2023-01-19 [2] Bioconductor
# sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.2.1)
# SingleCellExperiment * 1.20.1    2023-03-17 [2] Bioconductor
# sparseMatrixStats      1.10.0    2022-11-01 [2] Bioconductor
# statmod                1.5.0     2023-01-06 [2] CRAN (R 4.2.2)
# stringi                1.7.12    2023-01-11 [2] CRAN (R 4.2.2)
# stringr              * 1.5.0     2022-12-02 [2] CRAN (R 4.2.2)
# SummarizedExperiment * 1.28.0    2022-11-01 [2] Bioconductor
# tibble               * 3.2.1     2023-03-20 [2] CRAN (R 4.2.3)
# tidyr                * 1.3.0     2023-01-24 [2] CRAN (R 4.2.2)
# tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.2.1)
# tidyverse            * 2.0.0     2023-02-22 [2] CRAN (R 4.2.2)
# timechange             0.2.0     2023-01-11 [2] CRAN (R 4.2.2)
# tzdb                   0.3.0     2022-03-28 [2] CRAN (
