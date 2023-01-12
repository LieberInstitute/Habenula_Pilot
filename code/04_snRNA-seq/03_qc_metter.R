# January 12, 2023 - Bukola Ajanaku
# Applying QC metrics to filtered sce objects of the 7 Habenula samples.
# qrsh -l mem_free=50G,h_vmem=50G

# Loading relevant libraries
library("SingleCellExperiment")
library("jaffelab")
library("VariantAnnotation")
library("here")
library("ggplot2")
library("ggrepel")

# Loading filtered pre-QC sce object for all 7 Habenula samples
