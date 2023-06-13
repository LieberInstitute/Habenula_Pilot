## May 10, 2023 - Bukola Ajanaku
# Making stepwise plot using cleaned TSNE, cleaned num nuclei plots, and cleaning
# composition plots.
# qrsh -l mem_free=20G,h_vmem=20G

# loading relevant libraries
library(SummarizedExperiment)
library(here)
library(SingleCellExperiment)
library(DeconvoBuddies)
library(tidyverse)
library(tibble)
library(cowplot)
library(scater)

# loading sce object with dropped ambig cluster
load(file = here("processed-data", "99_paper_figs", "sce_objects", 
                 "sce_final_preHbdrop.RDATA"), verbose = TRUE)
table(sce$final_Annotations)

# dropping Excit.Neuron and OPC_noisy clusters
sce <- sce[, sce$final_Annotations != "OPC_noisy"]
sce <- sce[, sce$final_Annotations != "Excit.Neuron"]

# creating plot_dir
plot_dir <- here("plots", "99_paper_figs", "09_step_Wise_Plots")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

# sourcing official color palette 
source(file = here("code", "99_paper_figs", "source_colors.R"))
# bulk_colors and bulk_colors

############ PLOT 1: TSNE using sn annotations ##################################
# grabbing bulk annotations 
sce$bulkTypeSepHb <- sce$final_Annotations

# Combining into 5 glia, two thalamus, 1 broad LHb, and 1 broad MHb.
sce$bulkTypeSepHb[sce$bulkTypeSepHb %in% grep("^LHb\\.", unique(sce$bulkTypeSepHb), value = TRUE)] <-  'LHb'
sce$bulkTypeSepHb[sce$bulkTypeSepHb %in% grep("^MHb\\.", unique(sce$bulkTypeSepHb), value = TRUE)] <-  'MHb'

# check levels
table(sce$bulkTypeSepHb)

# cleaned TSNE with facet_wrap
TSNE <- plotReducedDim(sce, dimred = "TSNE") +
  geom_point(aes(color = sce$bulkTypeSepHb)) +
  scale_colour_manual(values = bulk_colors) +
  theme(legend.position = "none") +
  labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")

TSNE_facet <- plotReducedDim(sce, dimred = "TSNE") +
  geom_point(aes(color = sce$bulkTypeSepHb)) +
  scale_colour_manual(values = bulk_colors) +
  facet_wrap(~ sce$bulkTypeSepHb) +
  guides(color = guide_legend(title="Cell Type")) +
  labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")

pdf(file = here(plot_dir, "bulk_clean_TSNE.pdf"), width = 9, height = 5)
  plot_grid(TSNE, TSNE_facet)
dev.off()

# for One Drive
png(file = here(plot_dir, "mfigu_TSNE_by_CellType_CLEAN_FACETED_A.png"), width = 9, height = 5,
    units = "in", res = 1200)
  plot_grid(TSNE, TSNE_facet)
dev.off()

############ PLOT 2: TOTAL NUCLEI PLOT PER CT (bulk annotation) ################
# number of nuclei per cell type post drop
num_nuc <- as.data.frame(colData(sce)[,c("final_Annotations", 
                                         "bulkTypeSepHb", "Sample", "NeuN")]) |>
  group_by(Sample, bulkTypeSepHb, NeuN) |>
  mutate(n_nuc = n())

num_nuc_comp_plot <- num_nuc |>
  group_by(bulkTypeSepHb) |>
  summarize(tot_across_Samps = n()) |>
  ggplot(aes(x = bulkTypeSepHb, y = tot_across_Samps , fill = bulkTypeSepHb)) +
  geom_col() +
  geom_label(aes(label = tot_across_Samps),
             fill = "#FFFFFF",
             size = 4) +
  scale_fill_manual(values = bulk_colors) +
  theme_bw() +
  labs(y = "Number of Nuclei", fill = "Cell Type") +
  theme(axis.title.x = element_blank()) +
  theme(legend.position = "None", 
    axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())

pdf(file = here(plot_dir, "num_nuclei_post_clean.pdf"), width = 15, height = 6)
  num_nuc_comp_plot 
dev.off()

# for One Drive
png(file = here(plot_dir, "mfigu_num_Nuclei_by_CellType_CLEAN_B.png"), width = 7, height = 6,
    units = "in", res = 1200)
  num_nuc_comp_plot 
dev.off()


############ PLOT 3: TOTAL NUCLEI PLOT PER CT (bulk annotation) ################
prop_df <- as.data.frame(colData(sce)[, c("bulkTypeSepHb", "Sample")]) |>
  group_by(Sample, bulkTypeSepHb) |>
  summarize(n = n()) |>
  group_by(Sample) |>
  mutate(prop = n / sum(n))

comp_plot <- ggplot(prop_df, 
                        aes(x = Sample, y = prop, fill = bulkTypeSepHb)) +
  geom_col() +
  geom_text(
    aes(
    label = ifelse(prop > 0.02, format(round(prop, 3), 3), "")
    ), 
    size = 3,
    position = position_stack(vjust = 0.5),
    color = "white",
  ) +
  scale_fill_manual(values = bulk_colors) +
  theme_bw() +
  theme(legend.position = "None", 
        axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
  labs(y = "Proportion")


pdf(file = here(plot_dir, "comp_per_Sample_Bulk_Anno.pdf"), width = 10, height = 9)
  comp_plot
dev.off()

# for One Drive
png(file = here(plot_dir, "mfigu_sce_Sample_Comp_Plot_CLEAN_C.png"), width = 7, height = 6,
    units = "in", res = 1200)
  comp_plot
dev.off()

########## COMBINING PLOTS #####################################################
pdf(file = here(plot_dir, "stepWise_Plot_Post_Drop.pdf"), width = 10, height = 10)
plot_grid(
  TSNE, 
  TSNE_facet,
  num_nuc_comp_plot,
  comp_plot,
  labels = c("A", "", "B", "C")
)

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
# bluster                1.8.0     2022-11-01 [2] Bioconductor
# cli                    3.6.1     2023-03-23 [2] CRAN (R 4.2.3)
# cluster                2.1.4     2022-08-22 [3] CRAN (R 4.2.3)
# codetools              0.2-19    2023-02-01 [3] CRAN (R 4.2.3)
# colorout             * 1.2-2     2023-02-09 [1] Github (jalvesaq/colorout@79931fd)
# colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.2.2)
# cowplot              * 1.1.1     2020-12-30 [2] CRAN (R 4.2.1)
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
# ggbeeswarm             0.7.1     2022-12-16 [2] CRAN (R 4.2.2)
# ggplot2              * 3.4.2     2023-04-03 [2] CRAN (R 4.2.3)
# ggrepel                0.9.3     2023-02-03 [2] CRAN (R 4.2.2)
# glue                   1.6.2     2022-02-24 [2] CRAN (R 4.2.1)
# gridExtra              2.3       2017-09-09 [2] CRAN (R 4.2.1)
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
# scater               * 1.26.1    2022-11-13 [2] Bioconductor
# scran                  1.26.2    2023-01-19 [2] Bioconductor
# scuttle              * 1.8.4     2023-01-19 [2] Bioconductor
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
# vipor                  0.4.5     2017-03-22 [2] CRAN (R 4.2.1)
# viridis                0.6.2     2021-10-13 [2] CRAN (R 4.2.1)
# viridisLite            0.4.2     2023-05-02 [1] 