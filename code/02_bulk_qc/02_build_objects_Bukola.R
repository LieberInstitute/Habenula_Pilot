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

# Loading rse objects after brain swap #########################################

# gene
load(here("preprocessed_data", "count_data_bukola",  
          "rse_gene_Roche_Habenula_qcAndAnnotated_n69.Rdata"))
# exons
load(here("preprocessed_data", "count_data_bukola",  
          "rse_exon_Roche_Habenula_qcAndAnnotated_n69.Rdata"))
# transncripts
load(here("preprocessed_data", "count_data_bukola",  
          "rse_tx_Roche_Habenula_qcAndAnnotated_n69.Rdata"))
# junctions
load(here("preprocessed_data", "count_data_bukola",  
          "rse_jx_Roche_Habenula_qcAndAnnotated_n69.Rdata"))

# Data dimensions ##############################################################
dim(rse_gene)
# [1] 58037    69
dim(rse_exon)
# [1] 571623     69
dim(rse_jx)
# [1] 477623     69
dim(rse_tx)
# [1] 198093     69

# Verify data integrity ########################################################
## Checks each feature (row) to see if there are any NA sample values. 

which(rowSums(is.na(assay(rse_gene)) | assay(rse_gene) == "") > 0) # none
which(rowSums(is.na(assay(rse_exon)) | assay(rse_exon) == "") > 0) # none
which(rowSums(is.na(assay(rse_tx)) | assay(rse_tx) == "") > 0) # none
which(rowSums(is.na(assay(rse_jx)) | assay(rse_jx) == "") > 0) # none


# Identifying percentage of zeroes in sample counts across features #############
(length(which(assay(rse_gene) == 0)) * 100) / (nrow(rse_gene)*ncol(rse_gene))
# [1] 45.30061 

(length(which(assay(rse_exon) == 0)) * 100) / (nrow(rse_exon)*ncol(rse_exon))
# [1] 22.9808

(length(which(assay(rse_tx) == 0)) * 100) / (nrow(rse_tx)*ncol(rse_tx))
# [1] 33.21449

(length(which(assay(rse_jx) == 0)) * 100) / (nrow(rse_jx)*ncol(rse_jx))
# [1] 49.20676


# Normalizing read counts by transforming to counts per million (cpm) via edgeR
# gene
assays(rse_gene, withDimnames=FALSE)$logcounts = 
  edgeR::cpm(calcNormFactors(rse_gene, method = "TMM"), log = TRUE, prior.count = 0.5)

# exon
assays(rse_exon, withDimnames=FALSE)$logcounts = 
  edgeR::cpm(calcNormFactors(rse_exon, method = "TMM"), log = TRUE, prior.count = 0.5)

# transcript (these are already in transcripts per million rather than count,
# so we simply scale it to log2(TPM + 0.5)
assays(rse_tx, withDimnames=FALSE)$logcounts = log2(assays(rse_tx)$tpm + 0.5)

# junction 
assays(rse_jx, withDimnames=FALSE)$logcounts = 
  edgeR::cpm(calcNormFactors(rse_jx, method = "TMM"), log = TRUE, prior.count = 0.5)

## Focusing on rse_gene ########################################################
# Computing QC metrics for 



## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

