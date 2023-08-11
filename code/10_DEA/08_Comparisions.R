library("here")
library("data.table")
library("VennDiagram")
library("dplyr")
library("sessioninfo")

out_plot <- here("plots", "10_DEA", "08_Comparisions")
if (!dir.exists(out_plot)) dir.create(out_plot)
out_data <- here("processed-data", "10_DEA", "08_Comparisions")
if (!dir.exists(out_data)) dir.create(out_data)



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
