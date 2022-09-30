# September 27th, 2022
# RSE Gene Violin Plots - Bukola Ajanaku
# qrsh -l mem_free=24G,h_vmem=100G

library("SummarizedExperiment")
library("GenomicRanges")
library("ggplot2")
library("edgeR")
library ("sessioninfo") # trying to incorporate this package
library("scater")
# library("jaffelab") INSTALL

# Loads merged rse_gene from Geo, includes data from Leo.
load("/dcs04/lieber/lcolladotor/dbDev_LIBD001/RNAseq/4Bukola/rse_gene.merged.curated.n5780.rda")

# Changes variabble to "rse" for easier calling.
rse <- rse_gene
rm(rse_gene)

# Filters for only Protocol = RiboZeroGold, age, and sex.
filtRSE <- rse[, rse$Protocol == "RiboZeroGold" & rse$Sex == "M" &
                 (rse$Age >= 20.00 & rse$Age <= 69.00)]

# Computing logcounts.[ERROR]
assays(filtRSE, withDimnames = FALSE)$logcounts <-
  edgeR::cpm(calcNormFactors(filtRSE, method = "TMM"), log = TRUE, prior.count = 0.5)


# Rename filtRSE's rownames by Symbols for RowData
rownames(filtRSE) <- rowData(filtRSE)$Symbol

# (1) Loading Louise's code (save it as R code)
# source(pwd for file)

# Running code
p <- my_plotExpression(filtRSE, genes = c("GPR151"), ct = "Region")

ggsave(p, filename = "/users/bsimbiat/Habenula_Bulk/plots/70_GPR151/70_trialPlot2.png")


 pdf("/users/bsimbiat/Habenula_Bulk/plots/70_GPR151/70_trialPlot.pdf")
 p
 dev.off()
