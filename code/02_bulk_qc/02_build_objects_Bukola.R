# November 17, 2022
# 02_build_objects_Bukola.R - Building objects using relevant QC metrics for QC
# analysis as per  smokingMouse pipeline.
# qrsh -l mem_free=20G,h_vmem=20G

library(SummarizedExperiment)
library(recount)
library(edgeR)
library(jaffelab)
library(here)
# library(WGCNA) not available 
library(scater)
# library(biomartr) not available
library(sessioninfo)

# Loading rse objects after brain swap
load(here("preprocessed_data", "count_data_bukola",  
          "rse_gene_Roche_Habenula_qcAndAnnotated_n69.Rdata")) # gene info
rse = rse_gene

# Data dimensions
dim(rse)
# [1] 58037    69

# Verify data integrity
## checks each gene (row) to see if there are any NA values for each sample. 
which(rowSums(is.na(assay(rse)) | assay(rse) == "") > 0)
## named integer(0). [GOOD]

# Identifying percetage of zeroes in sample counts across genes
(length(which(assay(rse) == 0)) * 100) / (nrow(rse)*ncol(rse))
# [1] 45.30061 

# Normalizing read counts by transforming to counts per million (cpm) via edgeR
assays(rse_gene, withDimnames=FALSE)$logcounts = 
  edgeR::cpm(calcNormFactors(rse_gene, method = "TMM"), log = TRUE, prior.count = 0.5)






## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

