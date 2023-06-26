library("here")
library("SummarizedExperiment")
library("sessioninfo")



############################# Load rse gene object ############################

load(
    here(
        "processed-data",
        "rse_objects",
        "rse_gene_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata"
    ),
    verbose = TRUE
)
# Loading objects:
#   rse_gene

###############################################################################


######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################

