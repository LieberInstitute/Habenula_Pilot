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


# 