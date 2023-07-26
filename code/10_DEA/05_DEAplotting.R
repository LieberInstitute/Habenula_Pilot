library("here")
library("data.table")
library("dplyr")
library("gplots")
library("scater")
library("EnhancedVolcano")
library("sessioninfo")

out_plot <- here("plots", "10_DEA", "04_DEA")



################## Load significant results from DEA objects ##################

dea_res <- list.files(
    here(
        "processed-data", "10_DEA", "04_DEA"
    ),
    pattern = "DEA_All*",
    full.names = TRUE
)

dea_res <- lapply(dea_res, fread, data.table = FALSE)


###############################################################################

######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################

