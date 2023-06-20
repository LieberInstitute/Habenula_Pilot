library("here")
library("SummarizedExperiment")
library("ggplot2")
library("cowplot")
library("GGally")
library("sessioninfo")

output_path <- here("plots", "08_bulk_snpPC", "02_explore_snpPCs")

############################ Load rse gene objects ############################

load(
    here(
        "processed-data",
        "02_bulk_qc",
        "count_data_bukola",
        "rse_gene_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata"
    ),
    verbose = TRUE
)
# Loading objects:
#   rse_gene

lobstr::obj_size(rse_gene)
# 30.17 MB

class(rse_gene)
# [1] "RangedSummarizedExperiment"
# attr(,"package")
# [1] "SummarizedExperiment"

snpPCs <- read.table(
    here(
        "processed-data",
        "08_bulk_snpPC",
        "habenula_genotypes_n69.snpPCs.tab"
    ),
    header = TRUE
)

lobstr::obj_size(snpPCs)
# 11.80 kB

class(snpPCs)
# [1] "data.frame"

###############################################################################


######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
