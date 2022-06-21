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




sce.all.hb$prelimCluster <- factor(clusters.k20)
pdf(here("plots","09_snRNA-seq_re-processed", "prelimCluster_habenulan_n7.pdf"))
plotReducedDim(sce.all.hb, dimred="TSNE", colour_by="prelimCluster")
dev.off()

# Save
save(sce.all.hb, pc.choice.hb,
     file=here("processed-data","09_snRNA-seq_re-processed","03_clustering.Rda"))

# sgejobs::job_single('03_clustering', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript 03_clustering.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
