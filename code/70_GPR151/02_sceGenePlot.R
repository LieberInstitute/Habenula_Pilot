# September 30th, 2022
# SCE Gene Violin Plots - Bukola Ajanaku
# qrsh -l mem_free=24G,h_vmem=100G

library ("sessioninfo")
Sys.sleep(60)

library("SummarizedExperiment")
library("GenomicRanges")
library("ggplot2")
library("edgeR")
library("scater")
library("utils")
library("here")

# Loading the processed sce object made my Josh. [Double check with Louise]
load("/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/09_snRNA-seq_re-processed/07_annotation.Rda")

# Renaming sce.all.hb as sce for the purpose of this script.
sce <- sce.all.hb

# Loading my_plotExpression:
source(here("code", "70_GPR151", "rseViolin.R"))

# Reformatting sce for my_plotExpression:
rownames(sce) <- rowData(sce)$gene_name

m <- my_plotExpression(sce, "GPR151", assay = "logcounts", ct = "cellType")

ggsave(m, filename = here("code", "70_GPR151", "70_sceViolin.png"),
          height = 15, width = 18, units = "in")


print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()