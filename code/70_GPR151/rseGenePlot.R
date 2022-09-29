# September 27th, 2022
# RSE Gene Violin Plots - Bukola Ajanaku
# qrsh -l mem_free=24G,h_vmem=100G

library("SummarizedExperiment")
library("GenomicRanges")
library("ggplot2")
library("edgeR")
library ("sessioninfo") # trying to incorporate this package
library("scater")

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

# Subsets rse for GPR151.
subGPR <- filtRSE[rowData(filtRSE)$Symbol == "GPR151", ]

# Makes violin plots for logcounts of GPR151 by brain region.
# Create df (Sample ID, Dx, Sex, Logcounts)
# Using second method: scater::plotExpression:

# Can't be used for rse
# plotExpression(subGPR,  rownames(levels(colData(subGPR)$Region)))

# This is a very sloppy way to do this because any shift in the rows will render
# all the further analysis negligble. Currently banking on no shift. Hoping for
# guidance on a cleaner way to do this.
logcounts <- t(assays(subGPR)$logcounts)
region <- as.matrix(colData(subGPR)$Region)

df <- cbind(region, logcounts)
df <- as.data.frame(df)

trialPlot <- ggplot(df, aes(x = region, y = logcounts))

pdf("/users/bsimbiat/Habenula_Bulk/plots/70_GPR151/70_trialPlot")
trialPlot
dev.off()


# expression_long <- reshape2::melt(as.matrix(assays(subGPR)[[assays]][Dx, ]))
