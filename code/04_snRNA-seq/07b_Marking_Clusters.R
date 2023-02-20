## 2/16/23 - Bukola Ajanaku
# Labeling clusters based on gene expression plots.
# qrsh -l mem_free=20G,h_vmem=20G

library("SingleCellExperiment")
library("here")

# loading sce object with clustered harmonized data
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
          "sce_post_clustering_with_logcounts.Rdata"))