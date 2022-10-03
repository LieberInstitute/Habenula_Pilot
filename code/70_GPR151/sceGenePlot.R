# September 30th, 2022
# SCE Gene Violin Plots - Bukola Ajanaku
# qrsh -l mem_free=24G,h_vmem=100G

library("SummarizedExperiment")
library("GenomicRanges")
library("ggplot2")
library("edgeR")
library ("sessioninfo") # trying to incorporate this package
library("scater")
library("utils")

# Loading the processed sce object made my Josh. [Double check with Louise]
load("/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/09_snRNA-seq_re-processed/07_annotation.Rda")

# Renaming sce.all.hb as sce for the purpose of this script.
sce <- sce.all.hb

# Loading my_plotExpression:
source("/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/code/70_GPR151/rseViolin.R")

# Reformatting sce for my_plotExpression:
rownames(sce) <- rowData(sce)$gene_name

m <- my_plotExpression(sce, "GPR151", assay = "logcounts", ct = "cellType")

ggsave(m, filename = "/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/plots/70_GPR151/70_sceViolin.png",
          height = 15, width = 18, units = "in")
