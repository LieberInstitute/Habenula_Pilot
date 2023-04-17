## April 6, 2023 - Bukola Ajanaku
# TSNEs by my celltypes and custom color pallete 
# qrsh -l mem_free=20G,h_vmem=20G

library("SingleCellExperiment")
library("SummarizedExperiment")
library("here")
library("tidyverse")

# grabbing original sce (pre Hb drop but post OPC clean up)
load()