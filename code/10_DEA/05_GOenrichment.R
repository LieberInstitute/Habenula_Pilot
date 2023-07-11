library("here")
library("data.table")
library("dplyr")
library("clusterProfiler")
library("ggplot2")
library("sessioninfo")

out_plot <- here("plots", "10_DEA", "05_GOenrichment")
out_data <- here("processed-data", "10_DEA", "05_GOenrichment")



####################### Load tsv with all genes from DEA ######################

## For now I'm just gonna load these two models: qc-snpPCs-Hb and
## qc-totAssGene-snpPCs-Hb beacuse they had the best results

DE_qc_snpPCs_Hb <- fread(
    here(
        "processed-data",
        "10_DEA",
        "04_DEA",
        "DEA_AllGenes_qc-snpPCs-Hb.tsv"
    ),
    header = TRUE,
    sep = "\t",
    data.table = FALSE,
    stringsAsFactors = FALSE
)

DE_qc_totAssGene_snpPCs_Hb <- fread(
    here(
        "processed-data",
        "10_DEA",
        "04_DEA",
        "DEA_AllGenes_qc-totAssGene-snpPCs-Hb.tsv"
    ),
    header = TRUE,
    sep = "\t",
    data.table = FALSE,
    stringsAsFactors = FALSE
)

###############################################################################



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
