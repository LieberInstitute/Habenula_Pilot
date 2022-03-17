################################################################################
### LIBD pilot 10x snRNA-seq: DLPFC samples
###   **Region-specific analyses**
###   - R-batch job for detxn of optimal PC space with 'sce.all.hb' object
###         -> see '10x-pilot_region-specific_DLPFC_step02_clust-annot_MNTJan2020.R'
###            for setup of the SCE
### LAH 03May2021
################################################################################

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
#[1] 82

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
sample_prelimClusters[which(rowSums(sample_prelimClusters == 0) == 2),]



## check doublet score for each prelim clust
clusIndexes = splitit(sce.all.hb$prelimCluster)
prelimCluster.medianDoublet <- sapply(clusIndexes, function(ii){
  median(sce.all.hb$doubletScore[ii])
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
unname(clust.treeCut[order.dendrogram(dend)])
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
