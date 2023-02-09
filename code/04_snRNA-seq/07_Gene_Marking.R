## 2/9/23 - Bukola Ajanaku
# Annotating trails for my three top clustering methods.
# qrsh -l mem_free=20G,h_vmem=20G

library("dendextend")
library("dynamicTreeCut")
library("SingleCellExperiment")
library("batchelor")
library("scater")
library("scran")
library("uwot")
library("DropletUtils")
library("jaffelab")
library("Rtsne")
library("here")
library("utils")
library("sessioninfo")

# loading sce object with clustered harmonized data
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
     "sce_mid_clustering.Rdata"))

####### NOTE: ANNOTATING BY TWO TOP CLUSTERING METHODS #########################
# 1) "k_20_Erik" (Walk Trap 20: 22 Groups) Has Rand value (against regular ctErik) of 0.699.
# 2) "k_10_Erik" (WalkTrap 10: 36 Groups) Has Rand value (against regular ctErik) of 0.655.

