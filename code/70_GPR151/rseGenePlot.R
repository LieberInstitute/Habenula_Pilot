# September 27th, 2022
# RSE Gene Violin Plots - Bukola Ajanaku
# R

library("SummarizedExperiment")
library("GenomicRanges")
library("ggplot2")
library("edgeR")
library("jaffelab")

# Loads merged rse_gene from Geo, includes data from Leo.
load("/dcs04/lieber/lcolladotor/dbDev_LIBD001/RNAseq/4Bukola/rse_gene.merged.curated.n5780.rda")

# Changes variabble to "rse" for easier calling.
rse <- rse_gene
rm(rse_gene)

# Filters for only Protocol = RiboZeroGold (filters out BSP1 and some BSP2)
filtRSE <- rse[, rse$Protocol == "RiboZeroGold"]

# Computing logcounts.[ERROR]
assays(filtRSE, withDimnames = FALSE)$logcounts <-
    edgeR::cpm(calcNormFactors(filtRSE, method = "TMM"), log = TRUE, prior.count = 0.5)


# Subsets rse for GPR151.
subGPR <- filtRSE[rowData(filtRSE)$Symbol == "GPR151", ]


# Makes violin plots for logcounts of GPR151 by brain region.
