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


save(sce.all.hb, file = here("processed-data","09_snRNA-seq_re-processed","05_collapsedClustering.Rda"))

n_clusters <- length(levels(sce.all.hb$collapsedCluster))

## Add annotations, looking at marker gene expression
annotationTab.hab <- data.frame(collapsedCluster=c(1:n_clusters))
annotationTab.hab$cellType <- NA
annotationTab.hab$cellType[c(4)] <- paste0("Inhib")
annotationTab.hab$cellType[c(1,2,5,6)] <- paste0("Excit_", c("A","B","C","D"))
annotationTab.hab$cellType[c(9,10)] <- paste0("Astro_",c("A","B"))
annotationTab.hab$cellType[c(11)] <- c("Micro")
annotationTab.hab$cellType[c(6)] <- c("OPC")
annotationTab.hab$cellType[c(3)] <- c("Oligo")
