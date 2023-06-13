## May 3, 2023 - Bukola Ajanaku
# Running Bulk Deconvolution Plots for Top 10 gene markers! 
# qrsh -l mem_free=50G,h_vmem=50G

library(SummarizedExperiment)
library(here)
library(DeconvoBuddies)
library(SingleCellExperiment)
library(jaffelab)
library(dplyr)
library(BisqueRNA)
library(ggplot2)
library(tidyverse)

# loading final sce object 
load(file = here("processed-data", "99_paper_figs", "sce_objects", 
                 "official_final_sce.RDATA"), verbose = TRUE)
# sce

# loading final rse object (didn't change the location for this one)
load(here("processed-data", "02_bulk_qc", "count_data_bukola", 
          "rse_gene_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata"), verbose = TRUE)
# rse_gene

# creating plotting directory
plot_dir <- here("plots", "99_paper_figs", "10_bulk_Deconvo", "OPC_noisy")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

###### Adding bulk collapsed annotations to sce object #########################
sce$bulkTypeSepHb <- sce$final_Annotations
# making separated Hb (2)
sce$bulkTypeSepHb[sce$bulkTypeSepHb %in% grep("^LHb\\.", unique(sce$bulkTypeSepHb), value = TRUE)] <-  'LHb'
sce$bulkTypeSepHb[sce$bulkTypeSepHb %in% grep("^MHb\\.", unique(sce$bulkTypeSepHb), value = TRUE)] <-  'MHb'

# check!
table(sce$bulkTypeSepHb)
# Astrocyte   Endo Excit.Thal Inhib.Thal        LHb        MHb  Microglia 
# 538         38       1800       7612       2214        710        145 
# Oligo        OPC 
# 2178       1202 

###### Adding necessary symbols ################################################
sym_sce <- sce
rownames(sym_sce) <- rowData(sce)$Symbol

rownames(rse_gene) <- rowData(rse_gene)$Symbol

######## Pre-Bisque ############################################################
## remember, this is the broad analyses meaning that these annotations are solely
# for bulk deconvo

# Creating mean_ratios based on our specified annotations
ratios <- get_mean_ratio2(sym_sce,
                          cellType_col = "bulkTypeSepHb",
                          assay_name = "logcounts",
                          add_symbol = TRUE)

# Using the 1 vs All standard fold change for each gene x cell type
fc <- findMarkers_1vAll(sym_sce,
                        assay_name = "counts",
                        cellType_col = "bulkTypeSepHb",
                        add_symbol = FALSE,
                        mod = "~Sample",
                        verbose = TRUE
)

# combining the two to form marker_stats
marker_stats <- left_join(ratios, fc, by = c("gene", "cellType.target"))

# Random color scheme [NEED TO ESTABLISH MY OWN FOR THIS STEP]
cell_types <- levels(sym_sce$cellType)
# cell_colors <- create_cell_colors(cell_types = cell_types, pallet = "classic", split = "\\.", preview = TRUE)
# cell_colors errorr

# printing top 10 markers 
plot_marker_express_ALL(sym_sce,
                        marker_stats,
                        n_genes = 10,
                        rank_col = "rank_ratio",
                        anno_col = "anno_ratio",
                        cellType_col = "bulkTypeSepHb",
                        pdf_fn = here(plot_dir, "Top_10_Markers_OPC_noisy.pdf")
)

## copied directly from run_Bisque.R file in the bulk deconvo folder 
## creating marker_list of top 25 genes
marker_genes <- marker_stats |>
  filter(rank_ratio <= 25, gene %in% rownames(rse_gene)) |>
  pull(gene)

length(marker_genes)
# [1] 170

##### Running BISQUE ###########################################################
exp_set_bulk <- Biobase::ExpressionSet(assayData = assays(rse_gene[marker_genes,])$counts,
                                       phenoData=AnnotatedDataFrame(
                                         as.data.frame(colData(rse_gene))[c("BrNum")]))

exp_set_sce <- Biobase::ExpressionSet(assayData = as.matrix(assays(sym_sce[marker_genes,])$counts),
                                      phenoData=AnnotatedDataFrame(
                                        as.data.frame(colData(sym_sce))[,c("bulkTypeSepHb","Sample")]))

# checking for nuclei with 0 marker expression
zero_cell_filter <- colSums(exprs(exp_set_sce)) != 0
message("Exclude ", sum(!zero_cell_filter), " cells")
# Exclude 0 cells [yayyyyy!!!!]

est_prop <- ReferenceBasedDecomposition(bulk.eset = exp_set_bulk,
                                        sc.eset = exp_set_sce,
                                        cell.types = "bulkTypeSepHb",
                                        subject.names = "Sample",
                                        use.overlap = FALSE)

#### grabbed from my exploreBisque.R file in the bulk deconvo folder
# custom color scheme
color_bulk_clusters <- 
  c( "Oligo" = c("#A9A9A9"), # dark grey
     "OPC"= c("#7393B3"), # blue grey
     "Microglia" = c("#E5E4E2"), # platinum
     "Astrocyte" = c("#36454F"), # ash grey
     "Endo" = c("#848884"), # smoke
     "Inhib.Thal" = c('#2AAA8A'), # jungle green
     "Excit.Thal" = c("#478778"), # lincoln green
     "LHb" = c("#DE3163"), # cerise
     "MHb" = c("#00FFFF") # aqua
  )

# grabbing relevant phenotype info for bulk data
pd <- colData(rse_gene) |>
  as.data.frame() |>
  select(Sample = RNum, BrNum, PrimaryDx)

est_prop$bulk.props <- t(est_prop$bulk.props)
head(est_prop$bulk.props)

prop_long <- est_prop$bulk.props |>
  as.data.frame() |>
  rownames_to_column("Sample") |>
  tidyr::pivot_longer(!Sample, names_to = "cellType", values_to = "prop") |>
  left_join(pd) |>
  mutate(factor_CT = factor(cellType, levels = 
                              c("Astrocyte", "Endo", "Microglia", "Oligo", "OPC",
                                "Inhib.Thal", "Excit.Thal" , "MHb", "LHb")) ) |>
  arrange(factor_CT)

sum_Prop <- prop_long |>
  filter(cellType %in% c("LHb", "MHb")) |>
  group_by(BrNum) |>
  summarize(Hb_sum = sum(prop)) |>
  mutate(Br_Order = fct_reorder(BrNum, Hb_sum)) |>
  arrange(Br_Order)

prop_long <- left_join(prop_long, sum_Prop) |>
  arrange(Br_Order) 


## create composition bar plots
pdf(here(plot_dir, "bulk_Deconvo_Composition_OPC_noisy.pdf"), width = 21, height = 12)
plot_composition_bar(prop_long = prop_long, sample_col = "Br_Order",
                     x_col = "Br_Order", ct_col = "factor_CT") + 
  scale_fill_manual(values = color_bulk_clusters) +
  ggtitle("OPC_noisy")
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
# Biobase              * 2.58.0    2022-11-01 [2] Bioconductor
# BiocGenerics         * 0.44.0    2022-11-01 [2] Bioconductor
# BiocNeighbors          1.16.0    2022-11-01 [2] Bioconductor
# BiocParallel           1.32.6    2023-03-17 [2] Bioconductor
# BiocSingular           1.14.0    2022-11-01 [2] Bioconductor
# BisqueRNA            * 1.0.5     2021-05-23 [1] CRAN (R 4.2.3)
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
# fs                     1.6.1     2023-02-06 [2] CRAN (R 4.2.2)
# gargle                 1.3.0     2023-01-30 [2] CRAN (R 4.2.2)
# generics               0.1.3     2022-07-05 [2] CRAN (R 4.2.1)
# GenomeInfoDb         * 1.34.9    2023-02-02 [2] Bioconductor
# GenomeInfoDbData       1.2.9     2022-09-29 [2] Bioconductor
# GenomicRanges        * 1.50.2    2022-12-16 [2] Bioconductor
# ggplot2              * 3.4.2     2023-04-03 [2] CRAN (R 4.2.3)
# glue                   1.6.2     2022-02-24 [2] CRAN (R 4.2.1)
# googledrive            2.1.0     2023-03-22 [2] CRAN (R 4.2.3)
# gtable                 0.3.3     2023-03-21 [2] CRAN (R 4.2.3)
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.2.1)
# hms                    1.1.3     2023-03-21 [2] CRAN (R 4.2.3)
# igraph                 1.4.2     2023-04-07 [2] CRAN (R 4.2.3)
# IRanges              * 2.32.0    2022-11-01 [2] Bioconductor
# irlba                  2.3.5.1   2022-10-03 [2] CRAN (R 4.2.1)
# jaffelab             * 0.99.32   2023-02-09 [1] Github (LieberInstitute/jaffelab@7b7afe3)
# lattice                0.20-45   2021-09-22 [3] CRAN (R 4.2.3)
# lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.2.1)
# limma                  3.54.2    2023-02-28 [2] Bioconductor
# locfit                 1.5-9.7   2023-01-02 [2] CRAN (R 4.2.2)
# lubridate            * 1.9.2     2023-02-10 [2] CRAN (R 4.2.2)
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.2.1)
# MASS                   7.3-58.2  2023-01-23 [3] CRAN (R 4.2.3)
# Matrix                 1.5-4     2023-04-04 [2] CRAN (R 4.2.3)
# MatrixGenerics       * 1.10.0    2022-11-01 [2] Bioconductor
# matrixStats          * 1.0.0     2023-06-02 [1] CRAN (R 4.2.3)
# metapod                1.6.0     2022-11-01 [2] Bioconductor
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
# rlang                  1.1.1     2023-04-28 [1] CRAN (R 4.2.3)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.2.1)
# rsvd                   1.0.5     2021-04-16 [2] CRAN (R 4.2.1)
# S4Vectors            * 0.36.2    2023-02-26 [2] Bioconductor
# ScaledMatrix           1.6.0     2022-11-01 [2] Bioconductor
# scales                 1.2.1     2022-08-20 [2] CRAN (R 4.2.1)
# scran                  1.26.2    2023-01-19 [2] Bioconductor
# scuttle                1.8.4     2023-01-19 [2] Bioconductor
# segmented              1.6-4     2023-04-13 [1] CRAN (R 4.2.3)
# sessioninfo            1.2.2     2021-12-06 [2] CRAN (R 4.2.1)
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
# tzdb                   0.3.0     2022-03-28 [2] CRAN (R 4.2.1)
# utf8                   1.2.3     2023-01-31 [2] CRAN (R 4.2.2)
# vctrs                  0.6.2     2023-04-19 [1] CRAN (R 4.2.3)
# withr                  2.5.0     2022-03-03 [2] CRAN (R 4.2.1)
# XVector                0.38.0    2022-11-01 [2] Bioconductor
# zlibbioc               1.44.0    2022-11-01 [2] Bioconductor
# 
# [1] /users/bsimbiat/R/4.2.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/library

