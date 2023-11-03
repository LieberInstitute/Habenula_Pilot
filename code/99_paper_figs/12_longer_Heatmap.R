## June 14, 2023 - Bukola Ajanaku
# This is the exact code I used for the original progress report. Copy and pasted!
# Simply changed plotting directory!!!
# qrsh -l mem_free=30G,h_vmem=30G

library("SingleCellExperiment")
library("here")
library("sessioninfo")
library("jaffelab")
library("ggplot2")
library("tidyverse")
library("ComplexHeatmap")
library("spatialLIBD")
library("DeconvoBuddies")

# loading sce object
load(here("processed-data", "04_snRNA-seq",  "sce_objects", "sce_final.Rdata"),
     verbose = TRUE)
# sce_final 

# creating plot directory
plot_dir <- here("plots", "99_paper_figs", "12_longer_Heatmap")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

# sourcing official color palette 
source(file = here("code", "99_paper_figs", "source_colors.R"))
# bulk_colors and sn_colors


table(sce_final$final_Annotations)
# grabbing the resolution we need
sce_final <- sce_final[, !sce_final$final_Annotations %in% c("Astrocyte", "Endo", "Excit.Thal",
                                       "Inhib.Thal", "Microglia", "Oligo",
                                       "OPC") ]

# creating marker list
official_markers = list(
  "LHb.1_" = c("LINC02653","SYT2", "AC114316.1", "ONECUT2", "AC005906.2",
              "ATP8B1", "MMRN1", "FGF9", "GRIK4", "PRCP"), 
  "LHb.2_" = c("AC073071.1", "CRH", "AC012501.2", "LINC01572", "AC020611.2",
              "PTPRQ", 'KCNMB4', "MTMR8", "DTNBP1", "LRRN3"),
  "LHb.3_" = c ("ENTHD1","GRAP2", "MIR4500HG", "AC061958.1", "SYDE2", 
               "LGI2", "PCSK1", "MCOLN3", "AL121821.2", "LINC02112"),
  "LHb.4_" = c("TLE2", "CBLN1", "ANKRD34B", "SDC2", "GLRA1", "UBASH3B",
              "AC007656.2", "PCDH19", "KCNK2", "CLMP"),
  "LHb.5_" = c("LINC01619", "CCND2", "AC008574.1", "DIRAS3", "CDH7", "SEMA3D",
              "TRHDE", "KCTD1", "GNG12", "KCNC2"),
  "LHb.6_" = c("TACR3", "ESRP1", "AC011369.1", "KCNAB2", "AL354771.1", "CHRNA4",
              "MYO7A", "CHST8", "ANKRD50", "AP002989.1"),
  "LHb.7_" = c("AC008619.1", "LINC02296", "AC027312.1", "CLIC6", 'AC078845.1',
              "UBE2QL1", 'EAF2', "RSPO2", "AC009264.1","LINC01497"),
  "MHb.1_" = c("EXOC1L","F13A1", 'AC012499.1', "TAC1", "FREM2", "WDR72", 
              "TAC3", "PDE11A", 'AC107208.1', 'CCK'), 
  "MHb.2_" = c("CHAT", 'AC079760.2','LINC00616', "AC104170.1", "AC111198.1",
              "HRH3", "AC002069.2", 'SLC17A7', "PLCH2", 'NUP214'),
  "MHb.3_" = c("BHLHE22", "PKP2", "EBF3", "MIR548XHG", "FGF10-AS1", "FGF10",
              "FOXP4", "LINC01811", "AL136456.1", "PROX1-AS1")
)

# 'Hb' = c("GPR151", "CHRNB3", "POU4F1", "LINC01876", "CHRNA6")
# "MHb" = c("CHAT", "LINC01307", "NEUROD1", "CHRNB4", "LINC02143", "AC114321.1",
#           "AC104170.1", "AC079760.2", "AC024610.2", "AC022382.2"),
# "LHb" = c("HTR4", "BVES", "NRP1", "HTR2C", "CDH4", "COL25A1", "CBLN2", 
#           "CACNA1I", "WSCD2"),
# 'Neuron' = c('SYT1'),
# 'Exc_Neuron' = c('SLC17A6'), 
# 'Inh_Neuron' = c('GAD1')

#### check marker genes ####
marker_stats <- readxl::read_xlsx(here("plots", "99_paper_figs", "10c_snResolution_Top_Markers", "snResolution_top50MarkerGenes.xlsx"))
marker_stats[[1]] <- NULL

official_markers_tb <- tibble(cellType.target = gsub("_","",names(official_markers)), 
                              Symbol = official_markers) |>
  unnest(Symbol) 

marker_choice <- marker_stats |>
  right_join(official_markers_tb) |>
  select(cellType.target, Symbol, rank_ratio, rank_marker)
## all data driven 
marker_choice |> arrange(-rank_ratio)

#### Prep comple heatmap annotations ####
row_namers <- c( 
                "LHb.1",
                "LHb.2",
                "LHb.3",
                "LHb.4",
                "LHb.5",
                "LHb.6",
                "LHb.7",
                "MHb.1", 
                "MHb.2",
                "MHb.3"
)

# explicit color scheme 
color_official_markers = c(
  "LHb.1" = c("#0085af"),
  "LHb.2" = c("#0096FF"),
  "LHb.3" = c ("#89CFF0"), 
  "LHb.4" = c("#6F8FAF"), 
  "LHb.5" = c("#40E0D0"),
  "LHb.6" = c("#008080"), 
  "LHb.7" = c("#7DF9FF"), 
  "MHb.1" = c("#FF00FF"),
  "MHb.2" = c("#FAA0A0"), 
  "MHb.3" = c("#fa246a")
)

# "Hb" = c("#702963"), 
# "MHb" = c("#F33A6A"),
# "LHb" = c("#0000FF"), 
# 'Neuron' = c('#D8BFD8'),
# 'Exc_Neuron' = c("#9e4ad1"), 
# "Inh_Neuron" = c('#b5a2ff')

# checking colors
preview_colors <- function(cell_colors) {
  par(las = 2) # make label text perpendicular to axis
  par(mar = c(5, 8, 4, 2)) # increase y-axis margin.
  barplot(rep(1, length(cell_colors)),
          col = cell_colors,
          horiz = TRUE,
          axes = FALSE,
          names.arg = names(cell_colors)
  )
}

png(here(plot_dir, "gene_markers_color_pal.png"))
preview_colors(color_official_markers)
dev.off()

png(here(plot_dir, "clusters_color_pal.png"))
preview_colors(sn_colors)
dev.off()


# Pseudobulking to create compressed sce object
## faking the pseudobulking function out by setting sample as all same sample
sce_final$FakeSample <- "Br1011"
sce_final$RealSample <- sce_final$Sample
sce_final$Sample <- sce_final$FakeSample

set.seed(20220907) 
sce_simple_pb_snAnno3 <- registration_pseudobulk(sce_final, "final_Annotations", "Sample")


####### PLOTTING ###############################################################
# Plotting ComplexHeatmap
sce = sce_simple_pb_snAnno3
clusterMethod = "final_Annotations"
markerList = official_markers

# Replacing genes with symbols for heatmap (remember, this is pseudobulked data)
rownames(sce) <- rowData(sce)$Symbol

# renaming rownnames of colData(sce) based on row annotations
rownames(colData(sce)) <- paste(colData(sce)[, clusterMethod])

# reordering sce object for plottability
sce_reorder <- sce[unlist(markerList) , row_namers]

# Making data frame of genes we are interested in annd their general classification
markTable <- as.data.frame(unlist(markerList)) |> 
  rownames_to_column("cellType") |>
  rename(gene = `unlist(markerList)`) |>
  mutate(cellType = ss(cellType, "_", 1)) |>
  filter(gene %in% rowData(sce_reorder)$Symbol)

# getting z scores
marker_z_score <- scale(t(logcounts(sce_reorder)))
# corner(marker_z_score)

# heatmap columns annotation
column_ha <- HeatmapAnnotation(
  Gene_Marker = markTable$cellType,
  col = list(Gene_Marker = color_official_markers),
  show_legend = c(FALSE, FALSE)
  )

# grabbing the annotations per cluster from the sce_reorder object
clusterData <- as.data.frame(colData(sce_reorder)[,clusterMethod]) 
names(clusterData) <- "cellType"

# prepping the colors we want
# for cell type
# num_pal <- length(unique(clusterData$cellType))
# col_pal_ct <- grabColors(num_pal, start = 4)
# names(col_pal_ct) = unique(clusterData$cellType)

# heatmap row annotationn
row_ha <- rowAnnotation(
  Clusters = clusterData$cellType,
  col = list(Clusters = sn_colors),
  show_legend = FALSE
)

heatmapped <- Heatmap(marker_z_score,
                      col = circlize::colorRamp2(seq(-4, 4,  8/10),
                                                 rev(RColorBrewer::brewer.pal(11, "RdBu"))),
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      right_annotation = row_ha,
                      top_annotation = column_ha,
                      column_split = factor(markTable$cellType),
                      heatmap_legend_param = list(
                        title = c("Z_Score"),
                        border = "black"
                      ),
                      row_names_gp = grid::gpar(fontsize = 12),
                      column_names_gp = grid::gpar(fontsize = 12))

# printing 
pdf(here(plot_dir, "mainFigure_Hb_Markers_Heatmap.pdf"), width = 20, height = 6)
  heatmapped
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
# ! package                * version   date (UTC) lib source
# AnnotationDbi            1.60.2    2023-03-10 [2] Bioconductor
# AnnotationHub            3.6.0     2022-11-01 [2] Bioconductor
# attempt                  0.3.1     2020-05-03 [2] CRAN (R 4.2.1)
# beachmat                 2.14.2    2023-04-07 [2] Bioconductor
# beeswarm                 0.4.0     2021-06-01 [2] CRAN (R 4.2.1)
# benchmarkme              1.0.8     2022-06-12 [2] CRAN (R 4.2.1)
# benchmarkmeData          1.0.4     2020-04-23 [2] CRAN (R 4.2.1)
# Biobase                * 2.58.0    2022-11-01 [2] Bioconductor
# BiocFileCache            2.6.1     2023-02-17 [2] Bioconductor
# BiocGenerics           * 0.44.0    2022-11-01 [2] Bioconductor
# BiocIO                   1.8.0     2022-11-01 [2] Bioconductor
# BiocManager              1.30.20   2023-02-24 [2] CRAN (R 4.2.2)
# BiocNeighbors            1.16.0    2022-11-01 [2] Bioconductor
# BiocParallel             1.32.6    2023-03-17 [2] Bioconductor
# BiocSingular             1.14.0    2022-11-01 [2] Bioconductor
# BiocVersion              3.16.0    2022-04-26 [2] Bioconductor
# Biostrings               2.66.0    2022-11-01 [2] Bioconductor
# bit                      4.0.5     2022-11-15 [2] CRAN (R 4.2.2)
# bit64                    4.0.5     2020-08-30 [2] CRAN (R 4.2.1)
# bitops                   1.0-7     2021-04-24 [2] CRAN (R 4.2.1)
# blob                     1.2.4     2023-03-17 [2] CRAN (R 4.2.3)
# bluster                  1.8.0     2022-11-01 [2] Bioconductor
# bslib                    0.4.2     2022-12-16 [2] CRAN (R 4.2.2)
# cachem                   1.0.8     2023-05-01 [1] CRAN (R 4.2.3)
# circlize                 0.4.15    2022-05-10 [2] CRAN (R 4.2.1)
# cli                      3.6.1     2023-03-23 [2] CRAN (R 4.2.3)
# clue                     0.3-64    2023-01-31 [2] CRAN (R 4.2.2)
# cluster                  2.1.4     2022-08-22 [3] CRAN (R 4.2.3)
# codetools                0.2-19    2023-02-01 [3] CRAN (R 4.2.3)
# colorout               * 1.2-2     2023-02-09 [1] Github (jalvesaq/colorout@79931fd)
# colorspace               2.1-0     2023-01-23 [2] CRAN (R 4.2.2)
# ComplexHeatmap         * 2.14.0    2022-11-01 [2] Bioconductor
# config                   0.3.1     2020-12-17 [2] CRAN (R 4.2.1)
# cowplot                  1.1.1     2020-12-30 [2] CRAN (R 4.2.1)
# crayon                   1.5.2     2022-09-29 [2] CRAN (R 4.2.1)
# curl                     5.0.0     2023-01-12 [2] CRAN (R 4.2.2)
# data.table               1.14.8    2023-02-17 [2] CRAN (R 4.2.2)
# DBI                      1.1.3     2022-06-18 [2] CRAN (R 4.2.1)
# dbplyr                   2.3.2     2023-03-21 [2] CRAN (R 4.2.3)
# DeconvoBuddies         * 0.99.0    2023-03-10 [1] Github (lahuuki/DeconvoBuddies@d5774e3)
# DelayedArray             0.24.0    2022-11-01 [2] Bioconductor
# DelayedMatrixStats       1.20.0    2022-11-01 [2] Bioconductor
# digest                   0.6.31    2022-12-11 [2] CRAN (R 4.2.2)
# doParallel               1.0.17    2022-02-07 [2] CRAN (R 4.2.1)
# dotCall64                1.0-2     2022-10-03 [2] CRAN (R 4.2.1)
# dplyr                  * 1.1.1     2023-03-22 [2] CRAN (R 4.2.3)
# dqrng                    0.3.0     2021-05-01 [2] CRAN (R 4.2.1)
# DropletUtils             1.18.1    2022-11-22 [2] Bioconductor
# DT                       0.27      2023-01-17 [2] CRAN (R 4.2.2)
# edgeR                    3.40.2    2023-01-19 [2] Bioconductor
# ellipsis                 0.3.2     2021-04-29 [2] CRAN (R 4.2.1)
# ExperimentHub            2.6.0     2022-11-01 [2] Bioconductor
# fansi                    1.0.4     2023-01-22 [2] CRAN (R 4.2.2)
# fastmap                  1.1.1     2023-02-24 [2] CRAN (R 4.2.2)
# fields                   14.1      2022-08-12 [2] CRAN (R 4.2.1)
# filelock                 1.0.2     2018-10-05 [2] CRAN (R 4.2.1)
# forcats                * 1.0.0     2023-01-29 [2] CRAN (R 4.2.2)
# foreach                  1.5.2     2022-02-02 [2] CRAN (R 4.2.1)
# fs                       1.6.1     2023-02-06 [2] CRAN (R 4.2.2)
# gargle                   1.3.0     2023-01-30 [2] CRAN (R 4.2.2)
# generics                 0.1.3     2022-07-05 [2] CRAN (R 4.2.1)
# GenomeInfoDb           * 1.34.9    2023-02-02 [2] Bioconductor
# GenomeInfoDbData         1.2.9     2022-09-29 [2] Bioconductor
# GenomicAlignments        1.34.1    2023-03-09 [2] Bioconductor
# GenomicRanges          * 1.50.2    2022-12-16 [2] Bioconductor
# GetoptLong               1.0.5     2020-12-15 [2] CRAN (R 4.2.1)
# ggbeeswarm               0.7.1     2022-12-16 [2] CRAN (R 4.2.2)
# ggplot2                * 3.4.2     2023-04-03 [2] CRAN (R 4.2.3)
# ggrepel                  0.9.3     2023-02-03 [2] CRAN (R 4.2.2)
# GlobalOptions            0.1.2     2020-06-10 [2] CRAN (R 4.2.1)
# glue                     1.6.2     2022-02-24 [2] CRAN (R 4.2.1)
# golem                    0.4.0     2023-03-12 [2] CRAN (R 4.2.3)
# googledrive              2.1.0     2023-03-22 [2] CRAN (R 4.2.3)
# gridExtra                2.3       2017-09-09 [2] CRAN (R 4.2.1)
# gtable                   0.3.3     2023-03-21 [2] CRAN (R 4.2.3)
# HDF5Array                1.26.0    2022-11-01 [2] Bioconductor
# here                   * 1.0.1     2020-12-13 [2] CRAN (R 4.2.1)
# hms                      1.1.3     2023-03-21 [2] CRAN (R 4.2.3)
# htmltools                0.5.5     2023-03-23 [2] CRAN (R 4.2.3)
# htmlwidgets              1.6.2     2023-03-17 [2] CRAN (R 4.2.3)
# httpuv                   1.6.9     2023-02-14 [2] CRAN (R 4.2.2)
# httr                     1.4.5     2023-02-24 [2] CRAN (R 4.2.2)
# igraph                   1.4.2     2023-04-07 [2] CRAN (R 4.2.3)
# interactiveDisplayBase   1.36.0    2022-11-01 [2] Bioconductor
# IRanges                * 2.32.0    2022-11-01 [2] Bioconductor
# irlba                    2.3.5.1   2022-10-03 [2] CRAN (R 4.2.1)
# iterators                1.0.14    2022-02-05 [2] CRAN (R 4.2.1)
# jaffelab               * 0.99.32   2023-02-09 [1] Github (LieberInstitute/jaffelab@7b7afe3)
# jquerylib                0.1.4     2021-04-26 [2] CRAN (R 4.2.1)
# jsonlite                 1.8.5     2023-06-05 [1] CRAN (R 4.2.3)
# KEGGREST                 1.38.0    2022-11-01 [2] Bioconductor
# later                    1.3.0     2021-08-18 [2] CRAN (R 4.2.1)
# lattice                  0.20-45   2021-09-22 [3] CRAN (R 4.2.3)
# lazyeval                 0.2.2     2019-03-15 [2] CRAN (R 4.2.1)
# lifecycle                1.0.3     2022-10-07 [2] CRAN (R 4.2.1)
# limma                    3.54.2    2023-02-28 [2] Bioconductor
# locfit                   1.5-9.7   2023-01-02 [2] CRAN (R 4.2.2)
# lubridate              * 1.9.2     2023-02-10 [2] CRAN (R 4.2.2)
# magick                   2.7.4     2023-03-09 [2] CRAN (R 4.2.3)
# magrittr                 2.0.3     2022-03-30 [2] CRAN (R 4.2.1)
# maps                     3.4.1     2022-10-30 [2] CRAN (R 4.2.2)
# MASS                     7.3-58.2  2023-01-23 [3] CRAN (R 4.2.3)
# Matrix                   1.5-4     2023-04-04 [2] CRAN (R 4.2.3)
# MatrixGenerics         * 1.10.0    2022-11-01 [2] Bioconductor
# matrixStats            * 1.0.0     2023-06-02 [1] CRAN (R 4.2.3)
# memoise                  2.0.1     2021-11-26 [2] CRAN (R 4.2.1)
# metapod                  1.6.0     2022-11-01 [2] Bioconductor
# mime                     0.12      2021-09-28 [2] CRAN (R 4.2.1)
# munsell                  0.5.0     2018-06-12 [2] CRAN (R 4.2.1)
# nlme                     3.1-162   2023-01-31 [2] CRAN (R 4.2.2)
# paletteer                1.5.0     2022-10-19 [2] CRAN (R 4.2.1)
# pillar                   1.9.0     2023-03-22 [2] CRAN (R 4.2.3)
# pkgconfig                2.0.3     2019-09-22 [2] CRAN (R 4.2.1)
# plotly                   4.10.1    2022-11-07 [2] CRAN (R 4.2.2)
# png                      0.1-8     2022-11-29 [2] CRAN (R 4.2.2)
# promises                 1.2.0.1   2021-02-11 [2] CRAN (R 4.2.1)
# purrr                  * 1.0.1     2023-01-10 [2] CRAN (R 4.2.2)
# R.methodsS3              1.8.2     2022-06-13 [2] CRAN (R 4.2.1)
# R.oo                     1.25.0    2022-06-12 [2] CRAN (R 4.2.1)
# R.utils                  2.12.2    2022-11-11 [2] CRAN (R 4.2.2)
# R6                       2.5.1     2021-08-19 [2] CRAN (R 4.2.1)
# rafalib                * 1.0.0     2015-08-09 [1] CRAN (R 4.2.2)
# rappdirs                 0.3.3     2021-01-31 [2] CRAN (R 4.2.1)
# RColorBrewer             1.1-3     2022-04-03 [2] CRAN (R 4.2.1)
# Rcpp                     1.0.10    2023-01-22 [2] CRAN (R 4.2.2)
# RCurl                    1.98-1.12 2023-03-27 [2] CRAN (R 4.2.3)
# readr                  * 2.1.4     2023-02-10 [2] CRAN (R 4.2.2)
# rematch2                 2.1.2     2020-05-01 [2] CRAN (R 4.2.1)
# restfulr                 0.0.15    2022-06-16 [2] CRAN (R 4.2.1)
# rhdf5                    2.42.1    2023-04-07 [2] Bioconductor
# rhdf5filters             1.10.1    2023-03-24 [2] Bioconductor
# Rhdf5lib                 1.20.0    2022-11-01 [2] Bioconductor
# rjson                    0.2.21    2022-01-09 [2] CRAN (R 4.2.1)
# rlang                    1.1.1     2023-04-28 [1] CRAN (R 4.2.3)
# rprojroot                2.0.3     2022-04-02 [2] CRAN (R 4.2.1)
# Rsamtools                2.14.0    2022-11-01 [2] Bioconductor
# RSQLite                  2.3.1     2023-04-03 [2] CRAN (R 4.2.3)
# rsvd                     1.0.5     2021-04-16 [2] CRAN (R 4.2.1)
# rtracklayer              1.58.0    2022-11-01 [2] Bioconductor
# S4Vectors              * 0.36.2    2023-02-26 [2] Bioconductor
# R sass                     0.4.5     <NA>       [2] <NA>
#   ScaledMatrix             1.6.0     2022-11-01 [2] Bioconductor
# scales                   1.2.1     2022-08-20 [2] CRAN (R 4.2.1)
# scater                   1.26.1    2022-11-13 [2] Bioconductor
# scran                    1.26.2    2023-01-19 [2] Bioconductor
# scuttle                  1.8.4     2023-01-19 [2] Bioconductor
# segmented                1.6-4     2023-04-13 [1] CRAN (R 4.2.3)
# sessioninfo            * 1.2.2     2021-12-06 [2] CRAN (R 4.2.1)
# shape                    1.4.6     2021-05-19 [2] CRAN (R 4.2.1)
# shiny                    1.7.4     2022-12-15 [2] CRAN (R 4.2.2)
# shinyWidgets             0.7.6     2023-01-08 [2] CRAN (R 4.2.2)
# SingleCellExperiment   * 1.20.1    2023-03-17 [2] Bioconductor
# spam                     2.9-1     2022-08-07 [2] CRAN (R 4.2.1)
# sparseMatrixStats        1.10.0    2022-11-01 [2] Bioconductor
# SpatialExperiment      * 1.8.1     2023-03-05 [2] Bioconductor
# spatialLIBD            * 1.10.1    2022-12-01 [2] Bioconductor
# statmod                  1.5.0     2023-01-06 [2] CRAN (R 4.2.2)
# stringi                  1.7.12    2023-01-11 [2] CRAN (R 4.2.2)
# stringr                * 1.5.0     2022-12-02 [2] CRAN (R 4.2.2)
# SummarizedExperiment   * 1.28.0    2022-11-01 [2] Bioconductor
# tibble                 * 3.2.1     2023-03-20 [2] CRAN (R 4.2.3)
# tidyr                  * 1.3.0     2023-01-24 [2] CRAN (R 4.2.2)
# tidyselect               1.2.0     2022-10-10 [2] CRAN (R 4.2.1)
# tidyverse              * 2.0.0     2023-02-22 [2] CRAN (R 4.2.2)
# timechange               0.2.0     2023-01-11 [2] CRAN (R 4.2.2)
# tzdb                     0.3.0     2022-03-28 [2] CRAN (R 4.2.1)
# utf8                     1.2.3     2023-01-31 [2] CRAN (R 4.2.2)
# vctrs                    0.6.2     2023-04-19 [1] CRAN (R 4.2.3)
# vipor                    0.4.5     2017-03-22 [2] CRAN (R 4.2.1)
# viridis                  0.6.2     2021-10-13 [2] CRAN (R 4.2.1)
# viridisLite              0.4.2     2023-05-02 [1] CRAN (R 4.2.3)
# withr                    2.5.0     2022-03-03 [2] CRAN (R 4.2.1)
# XML                      3.99-0.14 2023-03-19 [2] CRAN (R 4.2.3)
# xtable                   1.8-4     2019-04-21 [2] CRAN (R 4.2.1)
# XVector                  0.38.0    2022-11-01 [2] Bioconductor
# yaml                     2.3.7     2023-01-23 [2] CRAN (R 4.2.2)
# zlibbioc                 1.44.0    2022-11-01 [2] Bioconductor
# 
# [1] /users/bsimbiat/R/4.2.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/library
# 
# R ── Package was removed from disk.


