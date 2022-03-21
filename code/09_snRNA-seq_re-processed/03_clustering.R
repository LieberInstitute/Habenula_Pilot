################################################################################
### LIBD pilot 10x snRNA-seq: DLPFC samples
###   **Region-specific analyses**
###   - R-batch job for detxn of optimal PC space with 'sce.all.hb' object
###         -> see '10x-pilot_region-specific_DLPFC_step02_clust-annot_MNTJan2020.R'
###            for setup of the SCE
### LAH 03May2021
################################################################################
library(dendextend)
library(dynamicTreeCut)
library(SingleCellExperiment)
library(scRNAseq)
library(batchelor)
library(EnsDb.Hsapiens.v86)
library(scater)
library(scran)
library(uwot)
library(DropletUtils)
library(jaffelab)
library(Rtsne)
library(here)
# ===


load(here("processed-data","09_snRNA-seq_re-processed","02_normalization.Rda"),
     verbose=TRUE)
    # sce.all.hb, chosen.hvgs.dlpfc

## PCA already done (interactively) - took top 100 PCs

## getClusteredPCs() to identify working PC space
pc.choice.hb <- getClusteredPCs(reducedDim(sce.all.hb))

# How many PCs should use in this space?
metadata(pc.choice.hb)$chosen
#[1] 57

## Plot n Clusters vs. d PCs
pdf(here("plots","09_snRNA-seq_re-processed", "PC_choice_habenulan_n7.pdf"))
plot(pc.choice.hb$n.pcs, pc.choice.hb$n.clusters,
     main=paste0("Combined Habnela (n=6) samples (d PCs choice = ", metadata(pc.choice.hb)$chosen, ")"))
abline(v=metadata(pc.choice.hb)$chosen, col="red", lty="dashed", lwd=0.8)
dev.off()


# Save
save(sce.all.hb, pc.choice.hb,
     file=here("processed-data","09_snRNA-seq_re-processed","03_clustering.Rda"))

# sgejobs::job_single('R-batchJob_DLPFC-n3_optimalPCselxn_LAH2021', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript R-batchJob_DLPFC-n3_optimalPCselxn_LAH2021.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()





set.seed(109)
sce.all.hb <- runTSNE(sce.all.hb, dimred="GLMPCA_approx")


## UMAP
set.seed(109)
sce.all.hb <- runUMAP(sce.all.hb, dimred="GLMPCA_approx")


# How do these look?
pdf(here("plots","09_snRNA-seq_re-processed", "ReducedDim_habenulan_n7.pdf"))
plotReducedDim(sce.all.hb, dimred="TSNE", colour_by="sample_short")
plotReducedDim(sce.all.hb, dimred="UMAP", colour_by="sample_short")
dev.off()

### Clustering: Two-step ======================================================
### Step 1: Perform graph-based clustering in this optimal PC space
#         - take k=20 NN to build graph
snn.gr <- buildSNNGraph(sce.all.hb, k=20, use.dimred="GLMPCA_approx")
clusters.k20 <- igraph::cluster_walktrap(snn.gr)$membership
table(clusters.k20)
    ## 107 prelim clusters



sce.all.hb$prelimCluster <- factor(clusters.k20)
pdf(here("plots","09_snRNA-seq_re-processed", "prelimCluster_habenulan_n7.pdf"))
plotReducedDim(sce.all.hb, dimred="TSNE", colour_by="prelimCluster")
dev.off()

# Is sample driving this 'high-res' clustering at this level?
(sample_prelimClusters <- table(sce.all.hb$prelimCluster, sce.all.hb$sample_short))  # (a little bit, but is typical)
  #  Br1092 Br1204 Br1469 Br1735 Br5555 Br5558 Br5639
  # 1      70      0    210    269     30     22    275
  # 2       3      0     22     19      7      6     33
  # 3    1054     24     39    158    487     29     47
  # 4       1      1   1053     22      1      2     19
  # 5       0      0    241      1      0      6      6
  # 6      15      0     98    158     17      5     87
  # 7       2      0     26    270      2      4     16
  # 8      13     85      0     27    154      0      7
  # 9     225     18      0      1    256      0      2
  # 10   1229     10      0   1880    937    354   1047
  # 11      0      0    232      0      0      0      0
  # 12      5      0     41    190      1      6     32
  # 13     41      3      0     26     56     22     14
  # 14      7      0    173    296      4     12    160
  # 15    164      0      0      4     29      0      0
  # 16     57      2      0    151    722      0      0
  # 17     20      0     84    113      2      2    167
  # 18      2    205      0      0      2      0      0
  # 19      3      0    188    132      0      1     19
  # 20     48      3      0     12    486      0      0
  # 21     17      0    115    143     31      5     97
  # 22      4      0     70    129      8      7     37
  # 23      0    207      0     10      3      0      0
  # 24     30      1      0     53    346      0      1
  # 25     11     29      1     52     14    426      8
  # 26      0     86      0      3      1      0      0
  # 27      1      0     69    375      4      4    320
  # 28      1     26      0      0      5      0      0
  # 29      0      9      0      0      1     26      0
  # 30      5      5      0      2     24      6      0
  # 31      0      0     78      2      0      0      6
  # 32    549      0      0    173    182      4   1001

sample_prelimClusters[which(rowSums(sample_prelimClusters == 0) == 2),]
  #    Br1092 Br1204 Br1469 Br1735 Br5555 Br5558 Br5639
  # 8      13     85      0     27    154      0      7
  # 9     225     18      0      1    256      0      2
  # 19      3      0    188    132      0      1     19
  # 24     30      1      0     53    346      0      1
  # 30      5      5      0      2     24      6      0
  # 32    549      0      0    173    182      4   1001


## check doublet score for each prelim clust

####FORMATED doublet data wrong. Will need to go back and generate later
clusIndexes = splitit(sce.all.hb$prelimCluster)
prelimCluster.medianDoublet <- sapply(clusIndexes, function(ii){
  median(sce.all.hb[[10]][ii])
}
)

summary(prelimCluster.medianDoublet)


hist(prelimCluster.medianDoublet)

## watch in clustering
prelimCluster.medianDoublet[prelimCluster.medianDoublet > 5]


#what are numbers here?
table(sce.all.hb$prelimCluster)[c(19, 32, 73)]
# 19 32 73
# 27 32  8



### Step 2: Hierarchical clustering of pseudo-bulked ("PB'd") counts with most robust normalization
  #         (as determined in: 'side-Rscript_testingStep2_HC-normalizn-approaches_wAmygData_MNTJan2020.R')
  #           ** That is, to pseudo-bulk (aka 'cluster-bulk') on raw counts, on all [non-zero] genes,
  #              normalize with `librarySizeFactors()`, log2-transform, then perform HC'ing

# Preliminary cluster index for pseudo-bulking
clusIndexes = splitit(sce.all.hb$prelimCluster)
prelimCluster.PBcounts <- sapply(clusIndexes, function(ii){
  rowSums(assays(sce.all.hb)$counts[ ,ii])
  }
)

    # And btw...
    table(rowSums(prelimCluster.PBcounts)==0)
# FALSE  TRUE
# 33824  2777

# Compute LSFs at this level
sizeFactors.PB.all  <- librarySizeFactors(prelimCluster.PBcounts)

# Normalize with these LSFs
geneExprs.temp <- t(apply(prelimCluster.PBcounts, 1, function(x) {log2(x/sizeFactors.PB.all + 1)}))

## Perform hierarchical clustering
dist.clusCollapsed <- dist(t(geneExprs.temp))
tree.clusCollapsed <- hclust(dist.clusCollapsed, "ward.D2")

dend <- as.dendrogram(tree.clusCollapsed, hang=0.2)



# Just for observation
pdf(here("plots","09_snRNA-seq_re-processed", "clust_dend_n7_habenula.pdf"))
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
#  0  1  2  3  4  5  6
# 10  6  4  4  3  3  2
unname(clust.treeCut[order.dendrogram(dend)])
# [1] 0 0 1 1 1 1 1 1 2 2 2 2 4 4 4 5 5 5 0 0 0 3 3 3 3 0 0 6 6 0 0 0
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
clust.treeCut[order.dendrogram(dend)][which(clust.treeCut[order.dendrogram(dend)]==0)] <- max(clust.treeCut)+c(1, 6, 2, 3, 4, 5, 5)

# 'Re-write', since there are missing numbers
# clust.treeCut[order.dendrogram(dend)] <- as.numeric(as.factor(clust2))
clust.treeCut[order.dendrogram(dend)] <- as.numeric(as.factor(clust.treeCut[order.dendrogram(dend)]))

## Define color pallet
cluster_colors <- unique(tableau20[clust.treeCut[order.dendrogram(dend)]])
names(cluster_colors) <- unique(clust.treeCut[order.dendrogram(dend)])
labels_colors(dend) <- cluster_colors[clust.treeCut[order.dendrogram(dend)]]

# Print for future reference
pdf(here("plots","09_snRNA-seq_re-processed", "perlimCluster-relationships_n7_habenula.pdf"), height = 9)
par(cex=0.6, font=2)
plot(dend, main="Habenula prelim-kNN-cluster relationships with collapsed assignments", horiz = TRUE)
abline(v = 325, lty = 2)
dev.off()


# Make reference for new cluster assignment
clusterRefTab.dlpfc <- data.frame(origClust=order.dendrogram(dend),
                                merged=clust.treeCut[order.dendrogram(dend)])


# Assign as 'collapsedCluster'
sce.all.hb$collapsedCluster <- factor(clusterRefTab.dlpfc$merged[match(sce.all.hb$prelimCluster, clusterRefTab.dlpfc$origClust)])
n_clusters <- length(levels(sce.all.hb$collapsedCluster))
# Print some visualizations:
pdf(here("plots","09_snRNA-seq_re-processed", "reducedDims_with-collapsedClusters_n7_habenula.pdf"))
plotReducedDim(sce.all.hb, dimred="PCA_corrected", ncomponents=5, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.all.hb, colour_by="sample_short", point_alpha=0.5)
plotTSNE(sce.all.hb, colour_by="cellType_Erik", point_alpha=0.5)
plotTSNE(sce.all.hb, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.all.hb, colour_by="doubletScore", point_alpha=0.5)
# And some more informative UMAPs
plotUMAP(sce.all.hb, colour_by="sample_short", point_alpha=0.5)
plotUMAP(sce.all.hb, colour_by="cellType_Erik", point_alpha=0.5)
plotUMAP(sce.all.hb, colour_by="collapsedCluster", point_alpha=0.5)
dev.off()

## Print marker genes for annotation
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/tableau_colors.rda", verbose = TRUE)

pdf("pdfs/revision/regionSpecific_DLPFC-n3_marker-logExprs_collapsedClusters_LAH2021.pdf", height=6, width=8)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpressionCustom(sce = sce.all.hb,
                         features = markers.mathys.custom[[i]],
                         features_name = names(markers.mathys.custom)[[i]],
                         anno_name = "collapsedCluster")
  )
}
dev.off()
# Assign as 'prelimCluster'
sce.all.hb$prelimCluster <- factor(clusters.k20)
plotReducedDim(sce.all.hb, dimred="TSNE", colour_by="prelimCluster")
