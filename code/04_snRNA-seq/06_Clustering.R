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
firstWay <- buildSNNGraph(sce, k = 10, use.dimred = "HARMONY")
lou <- igraph::cluster_louvain(firstWay)
sce$louvain <- paste0("Louvain", lou$membership)
table(sce$louvain)
    # Louvain1 Louvain10 Louvain11 Louvain12 Louvain13 Louvain14  Louvain2  Louvain3 
    # 1328       754       179      1274       542       893      2043      1485 
    # Louvain4  Louvain5  Louvain6  Louvain7  Louvain8  Louvain9 
    # 1436      2277       624      1538      1089      1620

##### Approach 2: WalkTrap Clustering (Erik) ###################################
g20 <- buildSNNGraph(sce, k = 20, use.dimred = 'HARMONY')
clust20 <- igraph::cluster_walktrap(g20)$membership
sce$k_20_Erik <- paste0("20wTrap", factor(clust20))
table(sce$k_20_Erik)
    # 20wTrap1 20wTrap10 20wTrap11 20wTrap12 20wTrap13 20wTrap14 20wTrap15 20wTrap16 
    # 1299       127       214      4747       296       192       151       347 
    # 20wTrap17 20wTrap18 20wTrap19  20wTrap2 20wTrap20 20wTrap21 20wTrap22  20wTrap3 
    # 211        34        75      1909        86       175        29      1867 
    # 20wTrap4  20wTrap5  20wTrap6  20wTrap7  20wTrap8  20wTrap9 
    # 569       450       359      2742       344       859 

g50 <- buildSNNGraph(sce, k=50, use.dimred = 'HARMONY')
clust50 <- igraph::cluster_walktrap(g50)$membership
colData(sce)$k_50_Erik <- paste0("50wTrap", factor(clust50))
table(colData(sce)$k_50_Erik)


### Saving data so far
save(sce, file = here("processed-data", "04_snRNA-seq", "sce_objects", 
                              "sce_mid_clustering.Rdata"))

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


