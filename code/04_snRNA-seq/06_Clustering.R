## 2/1/23 - Bukola Ajanaku
## Clustering post-harmony sce object (contains doublet scoring info) for clustering
## and future annotation.
## Based on:
# 1) https://github.com/LieberInstitute/Habenula_Bulk/blob/master/code/08_snRNA-seq_Erik/20210323_human_hb_neun.R
# 2) https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/code/03_build_sce/06_cluster.R
# 3) https://github.com/LieberInstitute/Habenula_Bulk/blob/master/code/09_snRNA-seq_re-processed/04_clustering.R
# 4) https://www.stephaniehicks.com/biocdemo/articles/Demo.html
# 5) OSCA workflows and multi-sample example sections: https://bioconductor.org/books/release/OSCA/book-contents.html
# qrsh -l mem_free=20G,h_vmem=20G

### TO DO:
# 2) Make plot by Sample and per Erik's cell type annotation
# 3) Create heatmap of comparing clusters

# Loading relevant libraries:
library("SingleCellExperiment")
library("jaffelab")
library("scater")
library("scran")
library("here")
library("sessioninfo")

# Loading post harmony sce object 
load(here("processed-data", "04_snRNA-seq", "sce_objects", "sce_harmony_by_Samp.Rdata"))
sce <- sce_corrbySamp
rm(sce_corrbySamp)

##### Approach 1: Louvain Clustering (Steph Hicks) #############################
set.seed(777)

firstWay <- buildSNNGraph(sce, k = 10, use.dimred = "HARMONY")
lou <- igraph::cluster_louvain(firstWay)
sce$louvain <- paste0("Louvain", lou$membership)
table(sce$louvain)
# Louvain1 Louvain10 Louvain11 Louvain12 Louvain13 Louvain14  Louvain2  Louvain3 
# 1328       756       178      1798       543       347      2254      1832 
# Louvain4  Louvain5  Louvain6  Louvain7  Louvain8  Louvain9 
# 1116      1855       620      1643      1089      1723


##### Approach 2: WalkTrap Clustering (Erik) ###################################
set.seed(1234)
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

set.seed(4321)
g50 <- buildSNNGraph(sce, k=50, use.dimred = 'HARMONY')
clust50 <- igraph::cluster_walktrap(g50)$membership
colData(sce)$k_50_Erik <- paste0("50wTrap", factor(clust50))
table(colData(sce)$k_50_Erik)
    # 50wTrap1 50wTrap10 50wTrap11 50wTrap12 50wTrap13 50wTrap14  50wTrap2  50wTrap3 
    # 110      4320       725       350       275       274       538       958 
    # 50wTrap4  50wTrap5  50wTrap6  50wTrap7  50wTrap8  50wTrap9 
    # 1826      2245       558       407      3174      1322 

set.seed(2565)
g10 <- buildSNNGraph(sce, k=10, use.dimred = 'HARMONY')
clust10 <- igraph::cluster_walktrap(g10)$membership
colData(sce)$k_10_Erik <-  paste0("10wTrap", factor(clust10))
table(colData(sce)$k_10_Erik)
    # 10wTrap1 10wTrap10 10wTrap11 10wTrap12 10wTrap13 10wTrap14 10wTrap15 10wTrap16 
    # 323        38      3218      1572      4254       216       149       499 
    # 10wTrap17 10wTrap18 10wTrap19  10wTrap2 10wTrap20 10wTrap21 10wTrap22 10wTrap23 
    # 202       169        97      1063       153       321        48       209 
    # 10wTrap24 10wTrap25 10wTrap26 10wTrap27 10wTrap28 10wTrap29  10wTrap3 10wTrap30 
    # 280        84       304        63        65        98      1326        65 
    # 10wTrap31 10wTrap32 10wTrap33 10wTrap34 10wTrap35 10wTrap36  10wTrap4  10wTrap5 
    # 22        37        39        17        38        17       138       179 
    # 10wTrap6  10wTrap7  10wTrap8  10wTrap9 
    # 686       337       371       385 

set.seed(14)
g5 <- buildSNNGraph(sce, k=5, use.dimred = 'HARMONY')
clust5 <- igraph::cluster_walktrap(g5)$membership
colData(sce)$k_5_Erik <- paste0("5wTrap", factor(clust5))
table(colData(sce)$k_5_Erik)
    # 5wTrap1  5wTrap10 5wTrap100 5wTrap101 5wTrap102 5wTrap103 5wTrap104 5wTrap105 
    # 61       275        20        10        14        23        14        11 
    # 5wTrap106 5wTrap107 5wTrap108  5wTrap11  5wTrap12  5wTrap13  5wTrap14  5wTrap15 
    # 7         7         8       337       413        97        70        48 
    # 5wTrap16  5wTrap17  5wTrap18  5wTrap19   5wTrap2  5wTrap20  5wTrap21  5wTrap22 
    # 94       142       127       133        47       112       120        33 
    # 5wTrap23  5wTrap24  5wTrap25  5wTrap26  5wTrap27  5wTrap28  5wTrap29   5wTrap3 
    # 120      1691       659        25       949        43        91        76 
    # 5wTrap30  ...


############# PLOTTING #########################################################
### Creating list to color by 
colorbyGroup <- c("louvain", "k_20_Erik", "k_50_Erik", "k_10_Erik", 
                  "ct_Erik")
## Getting rid of walkTrap 5 because we do not need 100+ grups 
# [1] "louvain"   "k_20_Erik" "k_50_Erik" "k_10_Erik"  

#### Previously included minibatch k means. No longer.


# Plotting harmonized UMAPS with colorbyGroup 
pdf(file = here("plots", "04_snRNA-seq", "06_Clustering", "UMAP_clustering_trails.pdf"),
    width = 7, height = 8)

lapply(colorbyGroup, function(m) {
#  plotUMAP(sce, colour_by = m) + theme(legend.position="bottom") + 
#    ggtitle(paste("Total Number of Groups =", length(table(colData(sce)[,m])))) 
    
})

dev.off()



# Plotting harmonized TSNE with colorbyGroup
pdf(file = here("plots", "04_snRNA-seq", "06_Clustering", "TSNE_clustering_trails2.pdf"),
    width = 7, height = 8)

#lapply(colorbyGroup, function(n) {
#  plotTSNE(sce, colour_by = n) + theme(legend.position="bottom") + 
#  ggtitle(paste("Total Number of Groups =", length(table(colData(sce)[,n]))))
  # # ggcells(sce, mapping = aes(x = TSNE.1, y = TSNE.2, colour = .data[[n]])) +
    ggcells(sce, mapping = aes(x = TSNE.1, y = TSNE.2, colour = "louvain")) +
        geom_point(size = 0.2, alpha = 0.3)  +
        my_theme +
        coord_equal() +
        labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2") +
        facet_wrap(~ Sample)
#})

dev.off()

### Saving sce object with clusters
save(sce, file = here("processed-data", "04_snRNA-seq", "sce_objects", 
                      "sce_mid_clustering.Rdata"))


## Issue with ^ sce save:
# table(duplicated(colnames(sce)))
# 
# FALSE  TRUE 
# 17048    34 

