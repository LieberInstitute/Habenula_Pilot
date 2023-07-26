library("here")
library("data.table")
library("dplyr")
library("gplots")
library("scater")
library("EnhancedVolcano")
library("sessioninfo")

out_plot <- here("plots", "10_DEA", "04_DEA")


######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################

