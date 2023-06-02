library("here")
library("SummarizedExperiment")
library("ggplot2")
library("scater")
library("rlang")
library("cowplot")
library("Hmisc")
library("lme4")
library("variancePartition")
library("reshape2")
library("pheatmap")
library("sessioninfo")



########################### Load rse filtered object ##########################

load(
    here(
        "processed-data", "10_DEA", "rse_gene_filt.Rdata"
    ),
    verbose = TRUE
)
# Loading objects:
#   rse_gene_filt

lobstr::obj_size(rse_gene_filt)
# 27.27 MB

class(rse_gene_filt)
# [1] "RangedSummarizedExperiment"
# attr(,"package")
# [1] "SummarizedExperiment"

###############################################################################



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################

