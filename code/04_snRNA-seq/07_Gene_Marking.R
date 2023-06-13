## 2/9/23 - Bukola Ajanaku
# Annotating trails for my three top clustering methods.
# qrsh -l mem_free=20G,h_vmem=20G

library("dendextend")
library("dynamicTreeCut")
library("SingleCellExperiment")
library("batchelor")
library("scater")
library("scran")
library("uwot")
library("DropletUtils")
library("jaffelab")
library("Rtsne")
library("here")
library("utils")
library("sessioninfo")
library("scuttle")

# loading sce object with clustered harmonized data
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
          "sce_post_clustering_with_logcounts.Rdata"))

####### NOTE: ANNOTATING BY TWO THREE CLUSTERING METHODS #######################
# Erik used Walktrap 50 and ended up with 14 groups.
# 1) Walktrap 50 (Rand: 0.661) has 14 groups.
# 2) Walktrap 10 (Rand: 0.654) has 37 groups.
# 3) Walktrap 20 (Rand: 0.632) has 23 groups.

####### SOURCING ###############################################################
# sourcing code from DLPFC Project (by Louise Huuki) 
# plots gene expression in a manner that renders images by fill rather than taking 
# up memory by plotting each point
  # generates expression plot 
source(here("code", "04_snRNA-seq", "sourcing", "custom_plotExpression.R"))

  # actually runs and plots markers 
source(here("code", "04_snRNA-seq", "sourcing", "my_plotMarkers.R"))

####### Marker gene list #######################################################
# [MUST HAVE AT LEAST TWO GENES PER MARKER FOR FUNCTION TO WORK CORRECTLY]
  # Erik's markers
    # markers.custom = list(
    #   'neuron' = c('SYT1', 'SNAP25'), #'SNAP25', 'GRIN1','MAP2'),
    #   'excitatory_neuron' = c('SLC17A6', 'SLC17A7'), # 'SLC17A8'),
    #   'inhibitory_neuron' = c('GAD1', 'GAD2'), #'SLC32A1'),
    #   'mediodorsal thalamus'= c('EPHA4','PDYN', 'LYPD6B', 'LYPD6', 'S1PR1', 'GBX2', 'RAMP3', 'COX6A2', 'SLITRK6', 'DGAT2'),
    #   'Hb neuron specific'= c('POU2F2','POU4F1','GPR151','CALB2'),#,'GPR151','POU4F1','STMN2','CALB2','NR4A2','VAV2','LPAR1'),
    #   'MHB neuron specific' = c('TAC1','CHAT','CHRNB4'),#'TAC3','SLC17A7'
    #   'LHB neuron specific' = c('HTR2C','MMRN1'),#'RFTN1'
    #   'oligodendrocyte' = c('MOBP', 'MBP'), # 'PLP1'),
    #   'oligodendrocyte_precursor' = c('PDGFRA', 'VCAN'), # 'CSPG4', 'GPR17'),
    #   'microglia' = c('C3', 'CSF1R'), #'C3'),
    #   'astrocyte' = c('GFAP', 'AQP4')
    # )

markers.custom = list(
  'neuron' = c('SYT1', 'SNAP25'), # 'GRIN1','MAP2'),
  'excitatory_neuron' = c('SLC17A6', 'SLC17A7'), # 'SLC17A8'),
  'inhibitory_neuron' = c('GAD1', 'GAD2'), #'SLC32A1'),
  'mediodorsal thalamus'= c('EPHA4','PDYN', 'LYPD6B', 'LYPD6', 'S1PR1', 'GBX2', 
                            'RAMP3', 'COX6A2', 'SLITRK6', 'DGAT2', "ADARB2"), # LYPD6B and ADARB2*, EPHA4 may not be best
  'Pf/PVT' = c("INHBA", "NPTXR"), 
  'Hb neuron specific'= c('POU4F1','GPR151'), #'STMN2','NR4A2','VAV2','LPAR1'), #CALB2 and POU2F2 were thrown out because not specific enough]
  'MHB neuron specific' = c('TAC1','CHAT','CHRNB4', "TAC3"),#'SLC17A7'
  'LHB neuron specific' = c('HTR2C','MMRN1', "ANO3"),#'RFTN1'
  'oligodendrocyte' = c('MOBP', 'MBP'), # 'PLP1'),
  'oligodendrocyte_precursor' = c('PDGFRA', 'VCAN'), # 'CSPG4', 'GPR17'),
  'microglia' = c('C3', 'CSF1R'), #'C3'),
  'astrocyte' = c('GFAP', 'AQP4'),
  "Endo/CP" = c("TTR", "FOLR1", "FLT1", "CLDN5")
)

#### PREPPING sce object to plot by gene expression ############################
# adding logcounts 
# sce <- logNormCounts(sce)
# 
# message("Start - saving data")
# Sys.time()
# 
# # saving before messing with row names for annotation purposes
# save(sce, file = here("processed-data", "04_snRNA-seq", "sce_objects",
#                         "sce_post_clustering_with_logcounts.Rdata"))
# message("Done - saving data")
# Sys.time()

# changing rownames for gene annotation purposes 
rownames(sce) <- rowData(sce)$Symbol


###### Plotting gene expression for walktrap method 10 (37 groups) #############
message("Start - Annotating wt10")
  Sys.time()
  
pdf10 <- here("plots", "04_snRNA-seq", "07_Gene_Marking", "wT_10_annotations_more_gran1.pdf")

my_plotMarkers(sce, marker_list = markers.custom, assay = "logcounts", 
               cat = "wT_10_Erik", fill_colors = NULL, pdf_fn = pdf10)

message("End - Annotating wt10")
Sys.time()

###### Plotting gene expression for walktrap method 20 (23 groups) #############
message("Start - Annotating wt20")
Sys.time()

pdf20 <- here("plots", "04_snRNA-seq", "07_Gene_Marking", "wT_20_annotations_more_gran1.pdf")

my_plotMarkers(sce, marker_list = markers.custom, assay = "logcounts", 
               cat = "wT_20_Erik", fill_colors = NULL, pdf_fn = pdf20)

message("End - Annotating wt20")
Sys.time()

###### Plotting gene expression for walktrap method 50 (14 groups) #############
message("Start - Annotating wt50")
  Sys.time()
  
  pdf50 <- here("plots", "04_snRNA-seq", "07_Gene_Marking", "wT_50_annotations_more_gran_1.pdf")
  
  my_plotMarkers(sce, marker_list = markers.custom, assay = "logcounts", 
                 cat = "wT_50_Erik", fill_colors = NULL, pdf_fn = pdf50)
  
  message("End - Annotating wt50")
Sys.time()



# sgejobs::job_single('07_Gene_Marking', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript 07_Gene_Marking.R")

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
# batchelor            * 1.14.1    2023-01-03 [2] Bioconductor
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
# DelayedArray           0.24.0    2022-11-01 [2] Bioconductor
# DelayedMatrixStats     1.20.0    2022-11-01 [2] Bioconductor
# dendextend           * 1.17.1    2023-03-25 [2] CRAN (R 4.2.3)
# dplyr                  1.1.1     2023-03-22 [2] CRAN (R 4.2.3)
# dqrng                  0.3.0     2021-05-01 [2] CRAN (R 4.2.1)
# DropletUtils         * 1.18.1    2022-11-22 [2] Bioconductor
# dynamicTreeCut       * 1.63-1    2016-03-11 [1] CRAN (R 4.2.2)
# edgeR                  3.40.2    2023-01-19 [2] Bioconductor
# fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.2.2)
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
# HDF5Array              1.26.0    2022-11-01 [2] Bioconductor
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.2.1)
# igraph                 1.4.2     2023-04-07 [2] CRAN (R 4.2.3)
# IRanges              * 2.32.0    2022-11-01 [2] Bioconductor
# irlba                  2.3.5.1   2022-10-03 [2] CRAN (R 4.2.1)
# jaffelab             * 0.99.32   2023-02-09 [1] Github (LieberInstitute/jaffelab@7b7afe3)
# lattice                0.20-45   2021-09-22 [3] CRAN (R 4.2.3)
# lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.2.1)
# limma                  3.54.2    2023-02-28 [2] Bioconductor
# locfit                 1.5-9.7   2023-01-02 [2] CRAN (R 4.2.2)
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.2.1)
# MASS                   7.3-58.2  2023-01-23 [3] CRAN (R 4.2.3)
# Matrix               * 1.5-4     2023-04-04 [2] CRAN (R 4.2.3)
# MatrixGenerics       * 1.10.0    2022-11-01 [2] Bioconductor
# matrixStats          * 1.0.0     2023-06-02 [1] CRAN (R 4.2.3)
# metapod                1.6.0     2022-11-01 [2] Bioconductor
# munsell                0.5.0     2018-06-12 [2] CRAN (R 4.2.1)
# nlme                   3.1-162   2023-01-31 [2] CRAN (R 4.2.2)
# pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.2.3)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.2.1)
# purrr                  1.0.1     2023-01-10 [2] CRAN (R 4.2.2)
# R.methodsS3            1.8.2     2022-06-13 [2] CRAN (R 4.2.1)
# R.oo                   1.25.0    2022-06-12 [2] CRAN (R 4.2.1)
# R.utils                2.12.2    2022-11-11 [2] CRAN (R 4.2.2)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.2.1)
# rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 4.2.2)
# RColorBrewer           1.1-3     2022-04-03 [2] CRAN (R 4.2.1)
# Rcpp                   1.0.10    2023-01-22 [2] CRAN (R 4.2.2)
# RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.2.3)
# ResidualMatrix         1.8.0     2022-11-01 [2] Bioconductor
# rhdf5                  2.42.1    2023-04-07 [2] Bioconductor
# rhdf5filters           1.10.1    2023-03-24 [2] Bioconductor
# Rhdf5lib               1.20.0    2022-11-01 [2] Bioconductor
# rlang                  1.1.1     2023-04-28 [1] CRAN (R 4.2.3)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.2.1)
# rsvd                   1.0.5     2021-04-16 [2] CRAN (R 4.2.1)
# Rtsne                * 0.16      2022-04-17 [2] CRAN (R 4.2.1)
# S4Vectors            * 0.36.2    2023-02-26 [2] Bioconductor
# ScaledMatrix           1.6.0     2022-11-01 [2] Bioconductor
# scales                 1.2.1     2022-08-20 [2] CRAN (R 4.2.1)
# scater               * 1.26.1    2022-11-13 [2] Bioconductor
# scran                * 1.26.2    2023-01-19 [2] Bioconductor
# scuttle              * 1.8.4     2023-01-19 [2] Bioconductor
# segmented              1.6-4     2023-04-13 [1] CRAN (R 4.2.3)
# sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.2.1)
# SingleCellExperiment * 1.20.1    2023-03-17 [2] Bioconductor
# sparseMatrixStats      1.10.0    2022-11-01 [2] Bioconductor
# statmod                1.5.0     2023-01-06 [2] CRAN (R 4.2.2)
# SummarizedExperiment * 1.28.0    2022-11-01 [2] Bioconductor
# tibble                 3.2.1     2023-03-20 [2] CRAN (R 4.2.3)
# tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.2.1)
# utf8                   1.2.3     2023-01-31 [2] CRAN (R 4.2.2)
# uwot                 * 0.1.14    2022-08-22 [2] CRAN (R 4.2.1)
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
# 
