## 2/1/23 - Bukola Ajanaku
## Clustering post-harmony sce object (contains doublet scoring info) for clustering
## and future annotation.
## Based on:
# 1) https://github.com/LieberInstitute/Habenula_Bulk/blob/master/code/08_snRNA-seq_Erik/20210323_human_hb_neun.R
# 2) https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/code/03_build_sce/06_cluster.R
# 3) https://github.com/LieberInstitute/Habenula_Bulk/blob/master/code/09_snRNA-seq_re-processed/04_clustering.R
# 4) https://www.stephaniehicks.com/biocdemo/articles/Demo.html
# 5) OSCA workflows and multi-sample example sections: https://bioconductor.org/books/release/OSCA/book-contents.html

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

##### Approach 1: Louvain Clustering (Steph Hicks) #############################
g <- buildSNNGraph(sce, k=10, use.dimred = "GLM-PCA")
lou <- igraph::cluster_louvain(g)
sce$louvain <- paste0("Louvain", lou$membership)
table(sce$louvain)

##### Approach 2: WalkTrap Clustering (Erik) ###################################
g20 <- buildSNNGraph(s3e.hb, k=20, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g20)$membership
colData(s3e.hb)$label <- factor(clust)
table(colLabels(s3e.hb))

g50 <- buildSNNGraph(s3e.hb, k=50, use.dimred = 'PCA')
clust50 <- igraph::cluster_walktrap(g50)$membership
colData(s3e.hb)$k_50_label <- factor(clust50)
table(colData(s3e.hb)$k_50_label)

g10 <- buildSNNGraph(s3e.hb, k=10, use.dimred = 'PCA')
clust10 <- igraph::cluster_walktrap(g10)$membership
colData(s3e.hb)$k_10_label <- factor(clust10)
table(colData(s3e.hb)$k_10_label)

g5 <- buildSNNGraph(s3e.hb, k=5, use.dimred = 'PCA')
clust5 <- igraph::cluster_walktrap(g5)$membership
colData(s3e.hb)$k_5_label <- factor(clust5)
table(colData(s3e.hb)$k_5_label)

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

##### Approach 4: DLPFC Clustering (Louise) ####################################
snn.gr <- buildSNNGraph(sce, k = 20, use.dimred = "HARMONY")
clusters <- igraph::cluster_walktrap(snn.gr)$membership

table(clusters)

