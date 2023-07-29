library("here")
library("dplyr")
library("purrr")
library("sessioninfo")

out_data <- here("processed-data", "12_GSEA")


######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################



