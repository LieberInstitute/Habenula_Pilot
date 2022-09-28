# September 27th, 2022
# RSE Gene Violin Plots - Bukola Ajanaku
# R

library("SummarizedExperiment")
library("GenomicRanges")
library("ggplot2")

# Loads "rse_gene" variable. Then change variabble to "rse" for easier calling.
load("/dcl02/lieber/ajaffe/Roche_Habenula/processed-data/01_bulk_speaqeasy/rse_gene_Roche_Habenula_qcAndAnnotated_n69.Rdata")
rse <- rse_gene
rm(rse_gene)


# Subsets rse for GPR151.
# Geo suggested a quicker method: grse <- rse[rowData(rse)$Symbol==('GPR151'), ]
findGPR <- rowRanges(rse)[rowRanges(rse)$Symbol == "GPR151"]
subGPR <- subsetByOverlaps(rse, findGPR)

# Filter out BSP1 and some BSP2 HPC (HIPPO) samples.


# Computes logcounts.


# Makes violin plots for logcounts of GPR151 by brain region.
