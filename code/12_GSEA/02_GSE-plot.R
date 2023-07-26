library("here")
library("sessioninfo")
library("dplyr")
library("purrr")
library("spatialLIBD")
library("ComplexHeatmap")
library("circlize")

out_plot <- here("plots", "12_GSEA")



######################### Reproducibility information #########################

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################

