################################################################################
### LIBD pilot 10x snRNA-seq: DLPFC samples
###   **Region-specific analyses**
###   - R-batch job for detxn of optimal PC space with 'sce.dlpfc' object
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
    # sce.dlpfc, chosen.hvgs.dlpfc

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
sce.all.hb <- runTSNE(sce.all.hb, dimred="PCA_opt")


## UMAP
set.seed(109)
sce.dlpfc <- runUMAP(sce.dlpfc, dimred="PCA_opt")


# How do these look?
plotReducedDim(sce.dlpfc, dimred="TSNE", colour_by="sampleID")
plotReducedDim(sce.dlpfc, dimred="UMAP", colour_by="sampleID")


### Clustering: Two-step ======================================================
### Step 1: Perform graph-based clustering in this optimal PC space
#         - take k=20 NN to build graph
snn.gr <- buildSNNGraph(sce.dlpfc, k=20, use.dimred="PCA_opt")
clusters.k20 <- igraph::cluster_walktrap(snn.gr)$membership
table(clusters.k20)
    ## 107 prelim clusters

# Is sample driving this 'high-res' clustering at this level?
(sample_prelimClusters <- table(sce.dlpfc$prelimCluster, sce.dlpfc$sampleID))  # (a little bit, but is typical)
sample_prelimClusters[which(rowSums(sample_prelimClusters == 0) == 2),]
# 39 - only 4 samples all from Br5207

# rbind the ref.sampleInfo[.rev]
ref.sampleInfo <- rbind(ref.sampleInfo, ref.sampleInfo.rev)

## check doublet score for each prelim clust
clusIndexes = splitit(sce.dlpfc$prelimCluster)
prelimCluster.medianDoublet <- sapply(clusIndexes, function(ii){
  median(sce.dlpfc$doubletScore[ii])
}
)

summary(prelimCluster.medianDoublet)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.01059  0.07083  0.14823  0.53264  0.30064 14.79144

hist(prelimCluster.medianDoublet)

## watch in clustering
prelimCluster.medianDoublet[prelimCluster.medianDoublet > 5]
# 19       32       73
# 14.79144  7.98099 10.20462

table(sce.dlpfc$prelimCluster)[c(19, 32, 73)]
# 19 32 73
# 27 32  8



### Step 2: Hierarchical clustering of pseudo-bulked ("PB'd") counts with most robust normalization
  #         (as determined in: 'side-Rscript_testingStep2_HC-normalizn-approaches_wAmygData_MNTJan2020.R')
  #           ** That is, to pseudo-bulk (aka 'cluster-bulk') on raw counts, on all [non-zero] genes,
  #              normalize with `librarySizeFactors()`, log2-transform, then perform HC'ing

# Preliminary cluster index for pseudo-bulking
clusIndexes = splitit(sce.dlpfc$prelimCluster)
prelimCluster.PBcounts <- sapply(clusIndexes, function(ii){
  rowSums(assays(sce.dlpfc)$counts[ ,ii])
  }
)

    # And btw...
    table(rowSums(prelimCluster.PBcounts)==0)
    # FALSE  TRUE
    # 29310  4228

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
par(cex=.6)
myplclust(tree.clusCollapsed, cex.main=2, cex.lab=1.5, cex=1.8)

dend %>%
  set("labels_cex", 0.8) %>%
  plot(horiz = TRUE)
abline(v = 325, lty = 2)

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
pdf("pdfs/revision/regionSpecific_DLPFC-n3_HC-prelimCluster-relationships_LAH2021.pdf", height = 9)
par(cex=0.6, font=2)
plot(dend, main="3x DLPFC prelim-kNN-cluster relationships with collapsed assignments", horiz = TRUE)
abline(v = 325, lty = 2)
dev.off()


# Make reference for new cluster assignment
clusterRefTab.dlpfc <- data.frame(origClust=order.dendrogram(dend),
                                merged=clust.treeCut[order.dendrogram(dend)])


# Assign as 'collapsedCluster'
sce.dlpfc$collapsedCluster <- factor(clusterRefTab.dlpfc$merged[match(sce.dlpfc$prelimCluster, clusterRefTab.dlpfc$origClust)])
n_clusters <- length(levels(sce.dlpfc$collapsedCluster))
# Print some visualizations:
pdf("pdfs/revision/regionSpecific_DLPFC-n3_reducedDims-with-collapsedClusters_LAH2021.pdf")
plotReducedDim(sce.dlpfc, dimred="PCA_corrected", ncomponents=5, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.dlpfc, colour_by="sampleID", point_alpha=0.5)
plotTSNE(sce.dlpfc, colour_by="protocol", point_alpha=0.5)
plotTSNE(sce.dlpfc, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.dlpfc, colour_by="sum", point_alpha=0.5)
plotTSNE(sce.dlpfc, colour_by="doubletScore", point_alpha=0.5)
# And some more informative UMAPs
plotUMAP(sce.dlpfc, colour_by="sampleID", point_alpha=0.5)
plotUMAP(sce.dlpfc, colour_by="collapsedCluster", point_alpha=0.5)
dev.off()

## Print marker genes for annotation
load(here("rdas","revision","markers.rda"), verbose = TRUE)

pdf("pdfs/revision/regionSpecific_DLPFC-n3_marker-logExprs_collapsedClusters_LAH2021.pdf", height=6, width=8)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpressionCustom(sce = sce.dlpfc,
                         features = markers.mathys.custom[[i]],
                         features_name = names(markers.mathys.custom)[[i]],
                         anno_name = "collapsedCluster")
  )
}
dev.off()
# Assign as 'prelimCluster'
sce.dlpfc$prelimCluster <- factor(clusters.k20)
plotReducedDim(sce.dlpfc, dimred="TSNE", colour_by="prelimCluster")
