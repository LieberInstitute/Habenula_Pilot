library(dendextend)
library(dynamicTreeCut)
library(SingleCellExperiment)
library(scRNAseq)
library(batchelor)
library(scater)
library(scran)
library(uwot)
library(DropletUtils)
library(jaffelab)
library(Rtsne)
library(here)
library(utils)
library(sessionInfo)
### Step 2: Hierarchical clustering of pseudo-bulked ("PB'd") counts with most robust normalization
  #         (as determined in: 'side-Rscript_testingStep2_HC-normalizn-approaches_wAmygData_MNTJan2020.R')
  #           ** That is, to pseudo-bulk (aka 'cluster-bulk') on raw counts, on all [non-zero] genes,
  #              normalize with `librarySizeFactors()`, log2-transform, then perform HC'ing
load(here("processed-data","09_snRNA-seq_re-processed","tableau_colors.rda"), verbose = TRUE)

load(here("processed-data","09_snRNA-seq_re-processed","03_clustering.Rda"))
# Preliminary cluster index for pseudo-bulking
clusIndexes = splitit(sce.all.hb$prelimCluster)
prelimCluster.PBcounts <- sapply(clusIndexes, function(ii){
  rowSums(assays(sce.all.hb)$counts[ ,ii])
  }
)

    # And btw...
    table(rowSums(prelimCluster.PBcounts)==0)
# FALSE  TRUE
# 33829  2772

# Compute LSFs at this level
sizeFactors.PB.all  <- librarySizeFactors(prelimCluster.PBcounts)

# Normalize with these LSFs
geneExprs.temp <- t(apply(prelimCluster.PBcounts, 1, function(x) {log2(x/sizeFactors.PB.all + 1)}))

## Perform hierarchical clustering
dist.clusCollapsed <- dist(t(geneExprs.temp))
tree.clusCollapsed <- hclust(dist.clusCollapsed, "ward.D2")

dend <- as.dendrogram(tree.clusCollapsed, hang=0.2)

# labels(dend)[grep(c("19|32|73"),labels(dend))] <- paste0(labels(dend)[grep(c("19|32|73"),labels(dend))], "*")

# Just for observation
pdf(here("plots","09_snRNA-seq_re-processed","dendrogram_hab_clusters.pdf"))
par(cex=.6)
myplclust(tree.clusCollapsed, cex.main=2, cex.lab=1.5, cex=1.8)

dend %>%
  set("labels_cex", 0.8) %>%
  plot(horiz = TRUE)
abline(v = 325, lty = 2)
dev.off()

clust.treeCut <- cutreeDynamic(tree.clusCollapsed, distM=as.matrix(dist.clusCollapsed),
                               minClusterSize=2, deepSplit=1, cutHeight=325)

table(clust.treeCut)
# clust.treeCut
# 0 1 2 3 4 5 6 7
# 4 5 4 4 3 2 2 2

unname(clust.treeCut[order.dendrogram(dend)])
#[1] 5 5 1 1 1 1 1 2 2 2 2 4 4 4 6 6 0 0 3 3 3 3 7 7 0 0
    ## Cutting at 250 looks good for the main neuronal branch, but a lot of glial
     #    prelim clusters are dropped off (0's)

    # # Cut at 400 for broad glia branch (will manually merge remaining dropped off)
    # glia.treeCut <- cutreeDynamic(tree.clusCollapsed, distM=as.matrix(dist.clusCollapsed),
    #                               minClusterSize=2, deepSplit=1, cutHeight=400)
    # unname(glia.treeCut[order.dendrogram(dend)])

    # Take those and re-assign to the first assignments

# clust <- clust.treeCut[order.dendrogram(dend)]
# clust2 <- name_zeros(clust, list(c(1,2), c(106,107)))
# unname(clust2)

# Add new labels to those prelimClusters cut off

# 'Re-write', since there are missing numbers
# clust.treeCut[order.dendrogram(dend)] <- as.numeric(as.factor(clust2))
clust.treeCut[order.dendrogram(dend)] <- as.numeric(as.factor(clust.treeCut[order.dendrogram(dend)]))

## Define color pallet
cluster_colors <- unique(tableau20[clust.treeCut[order.dendrogram(dend)]])
names(cluster_colors) <- unique(clust.treeCut[order.dendrogram(dend)])
labels_colors(dend) <- cluster_colors[clust.treeCut[order.dendrogram(dend)]]

# Print for future reference
pdf(here("plots","09_snRNA-seq_re-processed","regionSpecific_Hab-n7-prelimCluster-relationships.pdf"), height = 9)
par(cex=0.6, font=2)
plot(dend, main="7x Hab prelim-MNN-cluster relationships with collapsed assignments", horiz = TRUE)
abline(v = 325, lty = 2)
dev.off()


# Make reference for new cluster assignment
clusterRefTab.hab <- data.frame(origClust=order.dendrogram(dend),
                                merged=clust.treeCut[order.dendrogram(dend)])


# Assign as 'collapsedCluster'
sce.all.hb$collapsedCluster <- factor(clusterRefTab.hab$merged[match(sce.all.hb$prelimCluster, clusterRefTab.hab$origClust)])
n_clusters <- length(levels(sce.all.hb$collapsedCluster))
# Print some visualizations:
pdf(here("plots","09_snRNA-seq_re-processed","regionSpecific_HAB-n7_reducedDims-with-collapsedClusters.pdf"))
plotReducedDim(sce.all.hb, dimred="GLMPCA_MNN", ncomponents=5, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.all.hb, colour_by="sample_short", point_alpha=0.5)
plotTSNE(sce.all.hb, colour_by="high.mito", point_alpha=0.5)
plotTSNE(sce.all.hb, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.all.hb, colour_by="prelimCluster", point_alpha=0.5)
plotTSNE(sce.all.hb, colour_by="sum", point_alpha=0.5)
plotTSNE(sce.all.hb, colour_by="doubletScore", point_alpha=0.5)
plotTSNE(sce.all.hb, colour_by="cellType_Erik", point_alpha=0.5)

# And some more informative UMAPs
plotUMAP(sce.all.hb, colour_by="prelimCluster", point_alpha=0.5)
plotUMAP(sce.all.hb, colour_by="collapsedCluster", point_alpha=0.5)
plotUMAP(sce.all.hb, colour_by="sample_short", point_alpha=0.5)
plotUMAP(sce.all.hb, colour_by="high.mito", point_alpha=0.5)
plotUMAP(sce.all.hb, colour_by="sum", point_alpha=0.5)
plotUMAP(sce.all.hb, colour_by="doubletScore", point_alpha=0.5)
plotUMAP(sce.all.hb, colour_by="cellType_Erik", point_alpha=0.5)
dev.off()

save(sce.all.hb, file = here("processed-data","09_snRNA-seq_re-processed","05_collapsedClustering.Rda"))


print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
