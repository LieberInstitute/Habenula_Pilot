# November 17, 2022
# 02_build_objects_Bukola.R - Building objects using relevant QC metrics for QC
# analysis as per  smokingMouse pipeline.
# qrsh -l mem_free=20G,h_vmem=20G

library(SummarizedExperiment)
library(recount)
library(edgeR)
library(jaffelab)
library(here)
library(scater)
library(sessioninfo)
library(WGCNA) 
library(biomartr) 

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
# Computing QC metrics and adding back to column's metadata

# grabbing info on mitochondrial and ribosomal genes
subsets = list(Mito = which(seqnames(rse_gene)=="chrM"), 
             Ribo = grep("rRNA",rowData(rse_gene)$gene_type))

# addings qc stats based on the counts 
rse_gene <-addPerCellQC(rse_gene, subsets)

# grabbing relevant metadata 
pd = as.data.frame(colData(rse_gene))
gd = rowData(rse_gene)

# data filtration (getting rid of low expression values)
rse_gene_filt = rse_gene[which(filterByExpr(assay(rse_gene), 
      design = with(colData(rse_gene), model.matrix(~ AgeDeath + Flowcell + PrimaryDx)))),]

dim(rse_gene_filt)
# [1] 22756    69

# adding human gene symbols to replace MGI symbols 
rowData(rse_gene_filt)$MGI_Symbol<-rowData(rse_gene_filt)$Symbol

symbols<-biomart(genes  = rowData(rse_gene_filt)$ensemblID,
                 mart       = "ENSEMBL_MART_ENSEMBL",
                 dataset    = "mmusculus_gene_ensembl",
                 attributes = c("external_gene_name"),
                 filters    = "ensembl_gene_id")

# finding genes that did not have gene symbols
no_symbol = rowData(rse_gene_filt)$ensemblID[(! rowData(rse_gene_filt)$ensemblID 
            %in% symbols$ensembl_gene_id)]



# dropping irrelevant columns
drop = c("Brain.Region", "FQCbasicStats", "perBaseQual", "perTileQual",
         "GCcontent", "Ncontent", "SeqLengthDist", "SeqDuplication",
         "OverrepSeqs", "AdapterContent", "KmerContent", "SeqLength_R1",
         "perSeqQual", "perBaseContent", names(pd[,grepl("phred", names(pd))]),
         names(pd[,grepl("Adapter", names(pd))]), "SeqLength_R2", "bamFile",
         "trimmed", names(pd[,grepl("gene_", names(pd))]), "hasGenotype",
         "Age", "Race")
pd = pd[,!(names(pd)) %in% drop]


# Saving relevant variables for qc plotting
save(rse_gene, pd, gd, file = here("preprocessed_data", "count_data_bukola", 
                "built_objects_rse_gene_Roche_Habenula_n69.Rdata"))



## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

