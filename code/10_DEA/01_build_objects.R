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



################# Load deconvolution, SNP PCs and qSVa results ################

## Load deconvolution results
load(
    here(
        "processed-data",
        "06_deconvolution",
        "run_Bisque",
        "est_prop_split_Hb_annotations.RDATA"
    ),
    verbose = TRUE
)
# Loading objects:
#   est_prop


## Load SNP PCs
snpPCs <- read.table(
    here(
        "processed-data",
        "08_bulk_snpPC",
        "habenula_genotypes_filt.snpPCs.tab"
    ),
    header = TRUE
)

## Load qSVa

###############################################################################


######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################

