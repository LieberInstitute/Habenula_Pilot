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
# pd = as.data.frame(colData(rse_gene))
# gd = rowData(rse_gene)

# data filtration (getting rid of low expression values)
rse_gene_filt = rse_gene[which(filterByExpr(assay(rse_gene), 
      design = with(colData(rse_gene), model.matrix(~ AgeDeath + Flowcell + PrimaryDx)))),]

dim(rse_gene_filt)
# [1] 22756    69

# adding human gene symbols to replace MGI symbols 
# Note: gene symbols are abbrevation names for particular genes (established by 
# HUGO - Human Genome Organizatio)
rowData(rse_gene_filt)$MGI_Symbol<-rowData(rse_gene_filt)$Symbol

symbols<-biomart(genes  = rowData(rse_gene_filt)$ensemblID,
                 mart       = "ENSEMBL_MART_ENSEMBL",
                 dataset    = "mmusculus_gene_ensembl",
                 attributes = c("external_gene_name"),
                 filters    = "ensembl_gene_id")

# finding genes that did not have gene symbols
no_symbol = rowData(rse_gene_filt)$ensemblID[(! rowData(rse_gene_filt)$ensemblID 
            %in% symbols$ensembl_gene_id)]

# finding genes that have NA/empty symbols 
which_na_symbol = which(is.na(symbols$external_gene_name) | symbols$external_gene_name=="")
na_symbol <- symbols[which_na_symbol, 1]

# Compiling problematic gene IDs 
no_symbol = append(no_symbol, na_symbol)

# Removing problematic genes from symbols obj
symbols = symbols[-which_na_symbol,]

# Add Ensemble IDs for problematic genes
for (gene in no_symbol){
  
  MGI_symbol = rowData(rse_gene_filt)[which(rowData(rse_gene_filt)$ensemblID == gene), "MGI_Symbol"]
  if (! is.na(MGI_symbol)) {
    symbols[nrow(symbols) + 1,] = c(gene, MGI_symbol)
  }
  else {
    symbols[nrow(symbols)+1,]<-c(gene,gene)
  }
}

# Adding symbol to filtered rse_gene object while reserving original order of genes
symbols = symbols[match(rowData(rse_gene_filt)$ensemblID, symbols$ensembl_gene_id), ]
rse_gene = rse_gene_filt
rm(rse_gene_filt)
rowData(rse_gene)$Symbol = symbols$external_gene_name # external_gene_name is the MGI_symbol. 
save(rse_gene, file = here("processed-data", "02_bulk_qc", "count_data_bukola",  
                                "rse_gene_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata"))

## Repeating process for on non-gene rse objects ###############################
# doesn't require addPerCellQC as those QC metrics are only plotted for genotype

# EXONS
rse_exon_filt = rse_exon[which(filterByExpr(assay(rse_exon), 
                  design = with(colData(rse_exon), model.matrix(~ AgeDeath + Flowcell + PrimaryDx)))),]
dim(rse_exon_filt)
# [1] 356550     69

rowData(rse_exon_filt)$MGI_Symbol<-rowData(rse_exon_filt)$Symbol
symbols = biomart(genes  = rowData(rse_exon_filt)$ensemblID,
                 mart       = "ENSEMBL_MART_ENSEMBL",
                 dataset    = "mmusculus_gene_ensembl",
                 attributes = c("external_gene_name"),
                 filters    = "ensembl_gene_id")

# unique instead of colData (**)
no_symbol = unique(rowData(rse_exon_filt))$ensemblID[which(! unique(rowData(rse_exon_filt)$ensemblID) 
                                                  %in%  symbols$ensembl_gene_id)]
which_na_symbol = which(is.na(symbols$external_gene_name) | symbols$external_gene_name == "")
na_symbol = symbols[which_na_symbol, 1]
no_symbol = append(no_symbol, na_symbol)
symbols = symbols[-which_na_symbol,]

for (gene in no_symbol){
  
  MGI_symbol = unique(rowData(rse_exon_filt)[which(rowData(rse_exon_filt)$ensemblID==gene), "MGI_Symbol"])
  if (! (is.na(MGI_symbol) | length(MGI_symbol)==0)) {
    symbols[nrow(symbols)+1,]<-c(gene, MGI_symbol)
  }
  else {
    symbols[nrow(symbols)+1,]<-c(gene,gene)
  }
}

symbols = symbols[match(rowData(rse_exon_filt)$ensemblID, symbols$ensembl_gene_id), ]
rse_exon = rse_exon_filt
rm(rse_exon_filt)
rowData(rse_exon)$Symbol = symbols$external_gene_name # external_gene_name is the MGI_symbol. 
save(rse_exon, file = here("processed-data", "02_bulk_qc", "count_data_bukola",  
                                "rse_exon_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata"))

# JUNCTIONS
# filtration only. No symbols applicable since these are transcripts not genes
rse_jx_filt<-rse_jx[which(filterByExpr(assay(rse_jx), 
                    design=with(colData(rse_jx), model.matrix(~ AgeDeath + Flowcell + PrimaryDx)))),]
dim(rse_jx_filt)
# [1] 150926     69
save(rse_jx_filt, file = here("processed-data", "02_bulk_qc", "count_data_bukola",  
                                "rse_jx_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata"))

# TRANSCRIPTS
# Filtration involves creating potential cut offs since it was already in features per million
seeder = 12345432

pdf("processed-data/02_bulk_qc/count_data_bukola/tx_potential_cutoff.pdf")
  expression_cutoff(assays(rse_tx)$tpm, seed = seeder, k = 2)
dev.off()

# 2022-11-18 14:57:15 the suggested expression cutoff is 0.34
# percent_features_cut  samples_nonzero_cut 
# 0.42                 0.26 

cutoff = 0.34

rse_tx_filt = rse_tx[rowMeans(assays(rse_tx)$tpm) > cutoff,]
save(rse_tx_filt, file = here("processed-data", "02_bulk_qc", "count_data_bukola",  
                              "rse_tx_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata"))
dim(rse_tx_filt)
# [1] 82434    69

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()


