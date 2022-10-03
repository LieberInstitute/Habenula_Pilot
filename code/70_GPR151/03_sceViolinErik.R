# September 30th, 2022
# SCE Gene Violin Plots (using sce information from Erik, third plot for presentation) - Bukola Ajanaku
# qrsh -l mem_free=100G,h_vmem=100G

library("sessioninfo")
Sys.sleep(60)

library("SummarizedExperiment")
library("GenomicRanges")
library("ggplot2")
library("edgeR")
library("scater")
library("utils")
library("here")

# Loading sce data from Erik
load("/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/08_snRNA-seq_Erik/s3e_hb.rda")

# Renaming sce.all.hb as sce for the purpose of this script.
sce <- s3e.hb

# Loading my_plotExpression:
source(here("code", "70_GPR151", "rseViolin.R"))

# Reformatting sce for my_plotExpression:
rownames(sce) <- rowData(sce)$gene_name

# Running my_plotExpression:
m <- my_plotExpression(sce, genes = c("GPR151"), assay = "logcounts", ct = "cellType")

# Plotting
ggsave(m, filename = here("code", "70_GPR151", "70_sceViolinErik.png"),
       height = 15, width = 18, units = "in")

print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()