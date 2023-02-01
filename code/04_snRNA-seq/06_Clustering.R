## 2/1/23 - Bukola Ajanaku
## Clustering post-harmony sce object (contains doublet scoring info) for clustering
## and future annotation.
## Based on:
# 1) https://github.com/LieberInstitute/Habenula_Bulk/blob/master/code/08_snRNA-seq_Erik/20210323_human_hb_neun.R
# 2) https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/code/03_build_sce/06_cluster.R
# 3) https://github.com/LieberInstitute/Habenula_Bulk/blob/master/code/09_snRNA-seq_re-processed/04_clustering.R
# 4) https://www.stephaniehicks.com/biocdemo/articles/Demo.html
# 5) OSCA workflows and multi-sample example sections: https://bioconductor.org/books/release/OSCA/book-contents.html
# qrsh -l mem_free=50G,h_vmem=50G

# Loading relevant libraries:
library("SingleCellExperiment")
library("jaffelab")
library("scater")
library("scran")
library("here")
library("sessioninfo")
library("mbkmeans") # potentially not needed?

# Loading post harmony sce object 
load(here("processed-data", "04_snRNA-seq", "sce_objects", "sce_harmony_by_Samp.Rdata"))
sce <- sce_corrbySamp
rm(sce_corrbySamp)

##### Approach 1: Louvain Clustering (Steph Hicks) #############################
firstWay <- buildSNNGraph(sce, k=10, use.dimred = "HARMONY")
lou <- igraph::cluster_louvain(firstWay)
sce$louvain <- paste0("Louvain", lou$membership)
table(sce$louvain)

##### Approach 2: WalkTrap Clustering (Erik) ###################################
g20 <- buildSNNGraph(sce, k = 20, use.dimred = 'HARMONY')
clust20 <- igraph::cluster_walktrap(g20)$membership
colData(sce)$k_20_Erik <- factor(clust20)
table(colLabels(sce)$k_20_Erik)

g50 <- buildSNNGraph(sce, k=50, use.dimred = 'HARMONY')
clust50 <- igraph::cluster_walktrap(g50)$membership
colData(sce)$k_50_Erik <- factor(clust50)
table(colData(sce)$k_50_Erik)

g10 <- buildSNNGraph(sce, k=10, use.dimred = 'HARMONY')
clust10 <- igraph::cluster_walktrap(g10)$membership
colData(sce)$k_50_Erik <- factor(clust10)
table(colData(sce)$k_50_Erik)

g5 <- buildSNNGraph(sce, k=5, use.dimred = 'HARMONY')
clust5 <- igraph::cluster_walktrap(g5)$membership
colData(sce)$k_50_Erik <- factor(clust5)
table(colData(sce)$k_50_Erik)

##### Approach 3: Mini-batch k Means (Steph Hicks Example) #####################
k_list <- seq(5, 20)
km_res <- lapply(k_list, function(k) {
  mbkmeans(sce, clusters = k, 
           batch_size = 500,
           reduceMethod = "GLM-PCA",
           calc_wcss = TRUE)
})
wcss <- sapply(km_res, function(x) sum(x$WCSS_per_cluster))
plot(k_list, wcss, type = "b")


