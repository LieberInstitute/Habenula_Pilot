# September 27th, 2022
# RSE Gene Violin Plots - Bukola Ajanaku
# qrsh -l mem_free=100G,h_vmem=100G

library ("sessioninfo")
Sys.sleep(60)
library("SummarizedExperiment")
library("GenomicRanges")
library("ggplot2")
library("edgeR")
library("scater")
library("utils")
library("here")

# Loads merged rse_gene from Geo, includes data from Leo.
load("/dcs04/lieber/lcolladotor/dbDev_LIBD001/RNAseq/4Bukola/rse_gene.merged.curated.n5780.rda")

# Filters for only RiboZeroGold protocol, age, and sex.
filtRSE <- rse_gene[, rse_gene$Protocol == "RiboZeroGold" & rse_gene$Sex == "M" &
                      (rse_gene$Age >= 20.00 & rse_gene$Age <= 69.00)]

# Computing logcounts for filtRSE.
assays(filtRSE, withDimnames = FALSE)$logcounts <-
  edgeR::cpm(calcNormFactors(filtRSE, method = "TMM"), log = TRUE, prior.count = 0.5)

# Loading modified violin plot code for rse (my_plotExpression).
source(here("code", "70_GPR151", "rseViolin.R"))


# Prepping for my_plotExpression by renaming rownames by Symbols for RowData. Easy calling.
rownames(filtRSE) <- rowData(filtRSE)$Symbol

# Running code
p <- my_plotExpression(filtRSE, genes = c("GPR151"), ct = "Region")

ggsave(p, filename = here("code", "70_GPR151", "70_rseViolin.png"),
          height = 15, width = 18, units = "in")

print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()