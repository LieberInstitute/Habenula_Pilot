library("here")
library("SummarizedExperiment")
library("PCAtools")
library("stringr")
library("ggplot2")
library("cowplot")
library("GGally")
library("ComplexHeatmap")
library("circlize")
library("sessioninfo")

output_path <- here("plots", "10_DEA", "02_DataExploration")



############################# Load rse gene object ############################

load(
    here(
        "processed-data",
        "rse_objects",
        "rse_gene_Habenula_Pilot.rda"
    ),
)

lobstr::obj_size(rse_gene)
# 29.92 MB

class(rse_gene)
# [1] "RangedSummarizedExperiment"
# attr(,"package")
# [1] "SummarizedExperiment"

dim(rse_gene)
# [1] 22756    68

###############################################################################


######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
