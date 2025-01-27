## May 2, 2023 - Bukola Ajanaku
# Plotting cell-type expression pre and post drop per sample
# qrsh -l mem_free=50G,h_vmem=50G

# loading relevant libraries
library(SummarizedExperiment)
library(here)
library(SingleCellExperiment)
library(DeconvoBuddies)
library(tidyverse)
library(tibble)

# loading sce object with dropped ambig cluster
load(here("processed-data", "sce_objects", "sce_final_preHbdrop.RDATA"), verbose = TRUE)
##sce

table(sce$final_Annotations)
# Astrocyte         Endo Excit.Neuron   Excit.Thal   Inhib.Thal        LHb.1
# 538           38           51         1800         7612          201
# LHb.2        LHb.3        LHb.4        LHb.5        LHb.6        LHb.7
# 266          134          477           83           39         1014
# MHb.1        MHb.2        MHb.3    Microglia        Oligo          OPC
# 152          540           18          145         2178         1202
# OPC_noisy
# 594

# creating plot_dir
plot_dir <- here("plots", "99_paper_figs", "08_sce_Plot_Expression")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

# sourcing official color palette
source(file = here("code", "99_paper_figs", "source_colors.R"))
# bulk_colors and sn_colors


####### FINAL ANNOTATIONS LEVEL #################################################
#### get proportions before dropping ambig #####################################
# grabbing proportion information
prop_dirty_sn <- as.data.frame(colData(sce)[,
                                         c("final_Annotations", "Sample", "NeuN")]) |>
  group_by(Sample, final_Annotations, NeuN) |>
  summarize(n = n()) |>
  group_by(Sample) |>
  mutate(prop = n / sum(n))

# dropping the clusters we dropped
sce_drop <- sce[, sce$final_Annotations != "OPC_noisy"]
sce_drop <- sce_drop[, sce_drop$final_Annotations != "Excit.Neuron"]

#### proportions of nuclei using post-drop information #########################
pd_sn <- as.data.frame(colData(sce_drop))
table(pd_sn$final_Annotations)
# Astrocyte   Endo Excit.Thal Inhib.Thal      LHb.1      LHb.2      LHb.3
# 538         38       1800       7612        201        266        134
# LHb.4      LHb.5      LHb.6      LHb.7      MHb.1      MHb.2      MHb.3
# 477         83         39       1014        152        540         18
# Microglia Oligo        OPC
# 145       2178       1202

prop_clean_sn <- pd_sn[,c("final_Annotations", "Sample", "NeuN")] |>
  group_by(Sample, final_Annotations, NeuN) |>
  summarize(n = n()) |>
  group_by(Sample) |>
  mutate(prop = n / sum(n))

### combines prop_dirty and prop_clean
prop_ambig_plus_sn <- prop_dirty_sn |>
  mutate(Drop = "Pre-drop") |>
  bind_rows(prop_clean_sn |> mutate(Drop = "Post-drop"))

prop_ambig_plus_sn$Drop <- factor(prop_ambig_plus_sn$Drop,
                                  levels = c("Pre-drop", "Post-drop"))

# plots composition plot using prop_clean and prop_dirty
comp_plot_both_sn <- ggplot(data = prop_ambig_plus_sn, aes(x = Sample, y = prop,
                            fill = final_Annotations, group = Drop)) +
  geom_bar(stat = "identity") +
  geom_text(
    aes(
      label = ifelse(prop > 0.02, format(round(prop, 3), 3), "")
    ),
    size = 3,
    position = position_stack(vjust = 0.5),
    color = "black",
    family = "bold"
  ) +
  scale_fill_manual(values = c(sn_colors)) +
  labs(y = "Proportion", fill = "Cell Type") +
  facet_grid(Drop ~ NeuN, scales = "free", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color = "none", fill = guide_legend(ncol = 1,
                                             reverse = TRUE))
# pdf version
pdf(file = here(plot_dir, "sce_Comp_Plot_GRANULAR.pdf"), width = 7, height = 11)
  comp_plot_both_sn
dev.off()

# png version
png(file = here(plot_dir, "sce_Comp_Plot_GRANULAR.png"), width = 7, height = 11,
    units = "in", res = 1200)
  comp_plot_both_sn
dev.off()

### cell type sample breakdown - for reviews ##
prop_clean_sample <- pd_sn[,c("final_Annotations", "Sample", "NeuN")] |>
    group_by(Sample, final_Annotations, NeuN) |>
    summarize(n = n()) |>
    group_by(final_Annotations) |>
    mutate(prop = n / sum(n),
           final_Annotations = factor(final_Annotations, levels = names(sn_colors)))

comp_plot_sample <- ggplot(data = prop_clean_sample, aes(x = final_Annotations, y = prop,
                                                           fill = Sample)) +
    geom_bar(stat = "identity") +
    geom_text(
        aes(
            label = ifelse(prop > 0.02, format(round(prop, 2), 2), "")
        ),
        size = 3,
        position = position_stack(vjust = 0.5),
        color = "black"
    ) +
    labs(y = "Proportion", x = "Cell Type", fill = "Sample") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    guides(color = "none", fill = guide_legend(ncol = 1,
                                               reverse = TRUE))

ggsave(comp_plot_sample, file = here(plot_dir, "sce_Comp_Plot_Sample.png"), width = 7, height = 7)
ggsave(comp_plot_sample, file = here(plot_dir, "sce_Comp_Plot_Sample.pdf"), width = 7, height = 7)

####### BULK COLLAPSE LEVEL ####################################################
# creating bulk annotations level
sce$bulkTypeSepHb <- sce$final_Annotations
# Combining into 5 glia, two thalamus, 1 broad LHb, and 1 broad MHb.
sce$bulkTypeSepHb[sce$bulkTypeSepHb %in% grep("^LHb\\.", unique(sce$bulkTypeSepHb), value = TRUE)] <-  'LHb'
sce$bulkTypeSepHb[sce$bulkTypeSepHb %in% grep("^MHb\\.", unique(sce$bulkTypeSepHb), value = TRUE)] <-  'MHb'

# check levels
table(sce$bulkTypeSepHb)
    # Astrocyte         Endo Excit.Neuron   Excit.Thal   Inhib.Thal          LHb
    # 538           38           51         1800         7612         2214
    # MHb    Microglia        Oligo          OPC    OPC_noisy
    # 710          145         2178         1202          594

sce_drop$bulkTypeSepHb <- sce_drop$final_Annotations
sce_drop$bulkTypeSepHb[sce_drop$bulkTypeSepHb %in% grep("^LHb\\.", unique(sce_drop$bulkTypeSepHb), value = TRUE)] <-  'LHb'
sce_drop$bulkTypeSepHb[sce_drop$bulkTypeSepHb %in% grep("^MHb\\.", unique(sce_drop$bulkTypeSepHb), value = TRUE)] <-  'MHb'

#### get proportions before dropping ambig #####################################
prop_dirty_bulk <- as.data.frame(colData(sce)[,
                                c("bulkTypeSepHb", "Sample", "NeuN")]) |>
  group_by(Sample, bulkTypeSepHb, NeuN) |>
  summarize(n = n()) |>
  group_by(Sample) |>
  mutate(prop = n / sum(n))

#### proportions of nuclei using post-drop information #########################
# dirty
pd_bulk_dirty <- as.data.frame(colData(sce))
table(pd_bulk_dirty$bulkTypeSepHb)

prop_dirty_bulk <- pd_bulk_dirty[,c("bulkTypeSepHb", "Sample", "NeuN")] |>
  group_by(Sample, bulkTypeSepHb, NeuN) |>
  summarize(n = n()) |>
  group_by(Sample) |>
  mutate(prop = n / sum(n))

# clean
pd_bulk <- as.data.frame(colData(sce_drop))
table(pd_bulk$bulkTypeSepHb)

prop_clean_bulk <- pd_bulk[,c("bulkTypeSepHb", "Sample", "NeuN")] |>
  group_by(Sample, bulkTypeSepHb, NeuN) |>
  summarize(n = n()) |>
  group_by(Sample) |>
  mutate(prop = n / sum(n))

### combines prop_dirty and prop_clean
prop_ambig_plus_bulk <- prop_dirty_bulk |>
  mutate(Drop = "Pre-drop") |>
  bind_rows(prop_clean_bulk |> mutate(Drop = "Post-drop")) |>
  mutate(ct_levels = factor(bulkTypeSepHb, levels =
                              c("Excit.Neuron",
                                "Astrocyte",
                                "Endo",
                                "Microglia",
                                "Oligo",
                                "OPC_noisy",
                                "OPC",
                                "Inhib.Thal",
                                "Excit.Thal" ,
                                "MHb",
                                "LHb")) ) |>
  arrange(ct_levels)

prop_ambig_plus_bulk$Drop <- factor(prop_ambig_plus_bulk$Drop,
                                  levels = c("Pre-drop", "Post-drop"))


# plots composition plot using prop_clean and prop_dirty
comp_plot_both_bulk <- ggplot(data = prop_ambig_plus_bulk, aes(x = Sample,
                              y = prop, fill = bulkTypeSepHb, group = Drop)) +
  geom_bar(stat = "identity") +
  geom_text(
    aes(
      label = ifelse(prop > 0.02, format(round(prop, 3), 3), "")
    ),
    size = 3,
    position = position_stack(vjust = 0.5),
    color = "black"
  ) +
  scale_fill_manual(values = c(bulk_colors)) +
  labs(y = "Proportion", fill = "Cell Type") +
  facet_grid(Drop ~ NeuN, scales = "free", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color = "none", fill = guide_legend(ncol = 1,
                                             reverse = TRUE))
# pdf version
pdf(file = here(plot_dir, "sce_Comp_Plot_BROAD.pdf"), width = 7, height = 11)
  comp_plot_both_bulk
dev.off()

# png version
png(file = here(plot_dir, "sce_Comp_Plot_BROAD.png"), width = 7, height = 11,
    units = "in", res = 1200)
  comp_plot_both_bulk
dev.off()

# plots composition plot using prop_clean
comp_plot_clean <- ggplot(data = prop_clean_sn, aes(x = Sample,
                                                               y = prop, fill = final_Annotations)) +
  geom_bar(stat = "identity") +
  geom_text(
    aes(
      label = ifelse(prop > 0.02, format(round(prop, 3), 3), "")
    ),
    size = 3,
    position = position_stack(vjust = 0.5),
    color = "black"
  ) +
  scale_fill_manual(values = c(sn_colors)) +
  labs(y = "Proportion", fill = "Cell Type") +
  facet_wrap(~ NeuN, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  guides(color = "none", fill = guide_legend(ncol = 1,
                                             reverse = TRUE))
# pdf version
pdf(file = here(plot_dir, "sce_Comp_Plot_post-drop.pdf"), width = 3.5, height = 6.5)
comp_plot_clean
dev.off()

# png version
png(file = here(plot_dir, "sce_Comp_Plot_post-drop.png"), width = 3.5, height = 6.5,
    units = "in", res = 1200)
comp_plot_clean
dev.off()

# # plotting total nuclei information per sample
# barplot_n_nuc_bulk <- ggplot(prop_ambig_plus_bulk,
#   aes(x = Sample, y = n, fill = bulkTypeSepHb)) +
#   geom_col() +
#   geom_text(aes(label = n), size = 2.5) +
#   scale_fill_manual(values = bulk_colors) +
#   theme_bw() +
# #  theme(legend.position = "None", axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
#   labs(y = "Number of Nuclei") +
#   facet_grid(fct_rev(ambig) ~ NeuN, scales = "free", space = "free")
#
# pdf(file = here(plot_dir, "num_Nuc_Comp_Plot_bulkAnnoLEVEL.pdf"))
#   barplot_n_nuc_bulk
# dev.off()
#
# # plotting total nuclei information per sample
# sum_nuc_ambig_plus_prop <- prop_ambig_plus_bulk |>
#   group_by(ambig, Sample, bulkTypeSepHb) |>
#   summarize(n_across_samps = sum(n))
#
# barplot_n_nuc_bulk_tot <- ggplot(sum_nuc_ambig_plus_prop,
#          aes(x = bulkTypeSepHb, y = n_across_samps, fill = bulkTypeSepHb)) +
#   geom_col() +
#   geom_text(aes(label = n_across_samps), size = 2.5) +
#   scale_fill_manual(values = bulk_colors) +
#   theme_bw() +
#   #  theme(legend.position = "None", axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
#   labs(y = "Number of Nuclei") +
#   facet_wrap( ~ ambig, ncol = 1)
#
# pdf(file = here(plot_dir, "num_Nuc_Comp_Plot_bulkAnnoLEVEL_overall.pdf"), width = 10, height = 9)
#  barplot_n_nuc_bulk_tot
# dev.off()

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
# vctrs                  0.6.2     2023-


#
