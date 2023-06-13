# 2/15/23 - Bukola Ajanaku
# Reducing number of clusters rendered by walktrap 10 method.

library("here")
library("sessioninfo")
library("SingleCellExperiment")
library("dendextend")
library("dynamicTreeCut")
library("scater")
library("jaffelab")

# loading post clustered sce from 06_Clustering.R
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
          "sce_post_clustering.Rdata"))

message("Pseudobulk - ", Sys.time())
# using walktrap 10 as our prelim clusters
sce$prelimCluster <- sce$wT_10_Erik

# returns all indexes (column number in sce object) for each group 
clusIndexes <- splitit(sce$prelimCluster)

########### Check prelim clusters ##############################################

## Is sample driving this 'high-res' clustering at this level?
sample_prelimClusters <- table(sce$prelimCluster, sce$Sample)
sample_prelimClusters[which(rowSums(sample_prelimClusters != 0) == 1), ]
    # Br1092 Br1204 Br1469 Br1735 Br5555 Br5558 Br5639
    # 10wTrap_34      0      0     25      0      0      0      0
    # 10wTrap_37     17      0      0      0      0      0      0
  
#### check doublet score for each prelim clust ####
prelimCluster.medianDoublet <- sapply(clusIndexes, function(ii) {
  median(sce$doubletScore[ii])
})
summary(prelimCluster.medianDoublet)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.02999 0.25827 0.36485 0.43663 0.55054 1.03873 
# Nothing close to 5. Looks like no doublet driven clusters

########### PB ################################################################
# actual pseudo bulk. Takes counts for all indexes and makes it into one group
prelimCluster.PBcounts <- sapply(clusIndexes, function(ii) {
  rowSums(assays(sce)$counts[, ii])
})

# number of genes vs number of clusters
dim(prelimCluster.PBcounts)
  # [1] 33848    37

message("Get Lib Size Factors - ", Sys.time())
# Compute LSFs  at this level (adjusts size for each group)
sizeFactors.PB.all <- librarySizeFactors(prelimCluster.PBcounts)

message("Normalize - ", Sys.time())
# Normalize with these LSFs
geneExprs.temp <- t(apply(prelimCluster.PBcounts, 1, function(x) {
  log2(x / sizeFactors.PB.all + 1)
}))

message("Cluster Again - ", Sys.time())
## Perform hierarchical clustering
dist.clusCollapsed <- dist(t(geneExprs.temp))
tree.clusCollapsed <- hclust(dist.clusCollapsed, "ward.D2")

########### PLOTTING ###########################################################
plot_dir <- here("plots", "04_snRNA-seq", "08_Hierarchical_Clustering")
# dir.create(plot_dir)

# creating dendrogram (similar to treeplot)
dend <- as.dendrogram(tree.clusCollapsed, hang = 0.2)

# Print for future reference
pdf(here(plot_dir, "hClust_wTrap10_dendrogram.pdf"), height = 12, width = 12)
  par(cex = 0.6, font = 2)
  plot(dend, main = "hierarchical cluster dend", horiz = TRUE)
    # abline(v = 525, lty = 2)
dev.off()

# Cutting tree at chosen height
chosen_cut_height <- 150

clust.treeCut <- cutreeDynamic(tree.clusCollapsed,
                               distM = as.matrix(dist.clusCollapsed),
                               minClusterSize = 2, deepSplit = 1, cutHeight = chosen_cut_height
)

table(clust.treeCut)
unname(clust.treeCut[order.dendrogram(dend)])
unname(clust.treeCut[order.dendrogram(dend)]) == 0

## hmmm
    # > table(clust.treeCut)
    # clust.treeCut
    # 0  1  2 
    # 33  2  2 
    # > unname(clust.treeCut[order.dendrogram(dend)])
    # [1] 0 0 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    # > unname(clust.treeCut[order.dendrogram(dend)]) == 0
    # [1]  TRUE  TRUE FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    # [13]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE  TRUE
    # [25]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    # [37]  TRUE
    # > clust.treeCut
    # 2             2               3                             3           
    # 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 
    # > clust.treeCut[order.dendrogram(dend)]
    # 3 3                                   2 2                             
    # 0 0 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 

# Add new labels to those prelimClusters cut off
## just define as a max cluster for now
if (any(clust.treeCut[order.dendrogram(dend)] == 0)) {
  max_clust <- max(clust.treeCut)
  clust.treeCut[order.dendrogram(dend)][which(clust.treeCut[order.dendrogram(dend)] == 0)] <- max_clust + 1
  
  # 'Re-write', since there are missing numbers
  clust.treeCut[order.dendrogram(dend)] <- as.numeric(as.factor(clust.treeCut[order.dendrogram(dend)]))
  cluster_colors[as.character(max_clust + 1)] <- "black"
}

## Define HC color pallet
names(clust.treeCut) <- paste0("HC", numform::f_pad_zero(names(clust.treeCut)))
# cluster_colors <- DeconvoBuddies::create_cell_colors(cell_types = sort(unique(names(clust.treeCut))), pallet = "gg")

labels_colors(dend) <- cluster_colors[clust.treeCut[order.dendrogram(dend)]]

# Print for future reference
pdf(here(plot_dir, "dend_cut.pdf"), height = 11)
par(cex = 0.3, font = 1)
plot(dend, main = "DLPFC prelim-kNN-cluster relationships with collapsed assignments", horiz = TRUE)
abline(v = chosen_cut_height, lty = 2)
dev.off()


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
# dendextend           * 1.17.1    2023-03-25 [2] CRAN (R 4.2.3)
# dplyr                  1.1.1     2023-03-22 [2] CRAN (R 4.2.3)
# dynamicTreeCut       * 1.63-1    2016-03-11 [1] CRAN (R 4.2.2)
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
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.2.1)
# IRanges              * 2.32.0    2022-11-01 [2] Bioconductor
# irlba                  2.3.5.1   2022-10-03 [2] CRAN (R 4.2.1)
# jaffelab             * 0.99.32   2023-02-09 [1] Github (LieberInstitute/jaffelab@7b7afe3)
# lattice                0.20-45   2021-09-22 [3] CRAN (R 4.2.3)
# lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.2.1)
# limma                  3.54.2    2023-02-28 [2] Bioconductor
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.2.1)
# MASS                   7.3-58.2  2023-01-23 [3] CRAN (R 4.2.3)
# Matrix                 1.5-4     2023-04-04 [2] CRAN (R 4.2.3)
# MatrixGenerics       * 1.10.0    2022-11-01 [2] Bioconductor
# matrixStats          * 1.0.0     2023-06-02 [1] CRAN (R 4.2.3)
# munsell                0.5.0     2018-06-12 [2] CRAN (R 4.2.1)
# nlme                   3.1-162   2023-01-31 [2] CRAN (R 4.2.2)
# pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.2.3)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.2.1)
# purrr                  1.0.1     2023-01-10 [2] CRAN (R 4.2.2)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.2.1)
# rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 4.2.2)
# RColorBrewer           1.1-3     2022-04-03 [2] CRAN (R 4.2.1)
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
# segmented              1.6-4     2023-04-13 [1] CRAN (R 4.2.3)
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

