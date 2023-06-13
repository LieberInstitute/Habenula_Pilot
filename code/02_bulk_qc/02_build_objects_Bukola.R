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
colnames(rse_gene) <- rse_gene$RNum
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
colnames(rse_exon) <- rse_exon$RNum
save(rse_exon, file = here("processed-data", "02_bulk_qc", "count_data_bukola",  
                                "rse_exon_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata"))

# JUNCTIONS
# filtration only. No symbols applicable since these are transcripts not genes
rse_jx_filt<-rse_jx[which(filterByExpr(assay(rse_jx), 
                    design=with(colData(rse_jx), model.matrix(~ AgeDeath + Flowcell + PrimaryDx)))),]
dim(rse_jx_filt)
# [1] 150926     69
rse_jx_filt = rse_jx
rm(rse_jx_filt)
colnames(rse_jx) <- rse_jx$RNum
save(rse_jx, file = here("processed-data", "02_bulk_qc", "count_data_bukola",  
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
dim(rse_tx_filt)
# [1] 82434    69

rse_tx = rse_tx_filt
rm(rse_tx_filt)
colnames(rse_tx) <- rse_tx$RNum
save(rse_tx, file = here("processed-data", "02_bulk_qc", "count_data_bukola",  
                              "rse_tx_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata"))

## Reproducibility information
print('Reproducibility information:')
options(width = 120)
session_info()

# ─ Session info ──────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.2.3 Patched (2023-04-07 r84211)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2023-06-13
# pandoc   2.19.2 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/bin/ (via rmarkdown)
# 
# ─ Packages ──────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# AnnotationDbi          1.60.2    2023-03-10 [2] Bioconductor
# backports              1.4.1     2021-12-13 [2] CRAN (R 4.2.1)
# base64enc              0.1-3     2015-07-28 [2] CRAN (R 4.2.1)
# beachmat               2.14.2    2023-04-07 [2] Bioconductor
# beeswarm               0.4.0     2021-06-01 [2] CRAN (R 4.2.1)
# Biobase              * 2.58.0    2022-11-01 [2] Bioconductor
# BiocFileCache          2.6.1     2023-02-17 [2] Bioconductor
# BiocGenerics         * 0.44.0    2022-11-01 [2] Bioconductor
# BiocIO                 1.8.0     2022-11-01 [2] Bioconductor
# BiocNeighbors          1.16.0    2022-11-01 [2] Bioconductor
# BiocParallel           1.32.6    2023-03-17 [2] Bioconductor
# BiocSingular           1.14.0    2022-11-01 [2] Bioconductor
# biomaRt                2.54.1    2023-03-20 [2] Bioconductor
# biomartr             * 1.0.3     2023-05-07 [1] CRAN (R 4.2.3)
# Biostrings             2.66.0    2022-11-01 [2] Bioconductor
# bit                    4.0.5     2022-11-15 [2] CRAN (R 4.2.2)
# bit64                  4.0.5     2020-08-30 [2] CRAN (R 4.2.1)
# bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.2.1)
# blob                   1.2.4     2023-03-17 [2] CRAN (R 4.2.3)
# BSgenome               1.66.3    2023-02-16 [2] Bioconductor
# bumphunter             1.40.0    2022-11-01 [2] Bioconductor
# cachem                 1.0.8     2023-05-01 [1] CRAN (R 4.2.3)
# checkmate              2.1.0     2022-04-21 [2] CRAN (R 4.2.1)
# cli                    3.6.1     2023-03-23 [2] CRAN (R 4.2.3)
# cluster                2.1.4     2022-08-22 [3] CRAN (R 4.2.3)
# codetools              0.2-19    2023-02-01 [3] CRAN (R 4.2.3)
# colorout             * 1.2-2     2023-02-09 [1] Github (jalvesaq/colorout@79931fd)
# colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.2.2)
# crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.2.1)
# curl                   5.0.0     2023-01-12 [2] CRAN (R 4.2.2)
# data.table             1.14.8    2023-02-17 [2] CRAN (R 4.2.2)
# DBI                    1.1.3     2022-06-18 [2] CRAN (R 4.2.1)
# dbplyr                 2.3.2     2023-03-21 [2] CRAN (R 4.2.3)
# DelayedArray           0.24.0    2022-11-01 [2] Bioconductor
# DelayedMatrixStats     1.20.0    2022-11-01 [2] Bioconductor
# derfinder              1.32.0    2022-11-01 [2] Bioconductor
# derfinderHelper        1.32.0    2022-11-01 [2] Bioconductor
# digest                 0.6.31    2022-12-11 [2] CRAN (R 4.2.2)
# doParallel             1.0.17    2022-02-07 [2] CRAN (R 4.2.1)
# doRNG                  1.8.6     2023-01-16 [2] CRAN (R 4.2.2)
# downloader             0.4       2015-07-09 [2] CRAN (R 4.2.1)
# dplyr                  1.1.1     2023-03-22 [2] CRAN (R 4.2.3)
# dynamicTreeCut       * 1.63-1    2016-03-11 [1] CRAN (R 4.2.2)
# edgeR                * 3.40.2    2023-01-19 [2] Bioconductor
# evaluate               0.21      2023-05-05 [1] CRAN (R 4.2.3)
# fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.2.2)
# fastcluster          * 1.2.3     2021-05-24 [2] CRAN (R 4.2.1)
# fastmap                1.1.1     2023-02-24 [2] CRAN (R 4.2.2)
# filelock               1.0.2     2018-10-05 [2] CRAN (R 4.2.1)
# foreach                1.5.2     2022-02-02 [2] CRAN (R 4.2.1)
# foreign                0.8-84    2022-12-06 [3] CRAN (R 4.2.3)
# Formula                1.2-5     2023-02-24 [2] CRAN (R 4.2.2)
# fs                     1.6.1     2023-02-06 [2] CRAN (R 4.2.2)
# gargle                 1.3.0     2023-01-30 [2] CRAN (R 4.2.2)
# generics               0.1.3     2022-07-05 [2] CRAN (R 4.2.1)
# GenomeInfoDb         * 1.34.9    2023-02-02 [2] Bioconductor
# GenomeInfoDbData       1.2.9     2022-09-29 [2] Bioconductor
# GenomicAlignments      1.34.1    2023-03-09 [2] Bioconductor
# GenomicFeatures        1.50.4    2023-01-24 [2] Bioconductor
# GenomicFiles           1.34.0    2022-11-01 [2] Bioconductor
# GenomicRanges        * 1.50.2    2022-12-16 [2] Bioconductor
# GEOquery               2.66.0    2022-11-01 [2] Bioconductor
# ggbeeswarm             0.7.1     2022-12-16 [2] CRAN (R 4.2.2)
# ggplot2              * 3.4.2     2023-04-03 [2] CRAN (R 4.2.3)
# ggrepel                0.9.3     2023-02-03 [2] CRAN (R 4.2.2)
# glue                   1.6.2     2022-02-24 [2] CRAN (R 4.2.1)
# GO.db                  3.16.0    2022-09-28 [2] Bioconductor
# googledrive            2.1.0     2023-03-22 [2] CRAN (R 4.2.3)
# gridExtra              2.3       2017-09-09 [2] CRAN (R 4.2.1)
# gtable                 0.3.3     2023-03-21 [2] CRAN (R 4.2.3)
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.2.1)
# Hmisc                  5.0-1     2023-03-08 [2] CRAN (R 4.2.3)
# hms                    1.1.3     2023-03-21 [2] CRAN (R 4.2.3)
# htmlTable              2.4.1     2022-07-07 [2] CRAN (R 4.2.1)
# htmltools              0.5.5     2023-03-23 [2] CRAN (R 4.2.3)
# htmlwidgets            1.6.2     2023-03-17 [2] CRAN (R 4.2.3)
# httr                   1.4.5     2023-02-24 [2] CRAN (R 4.2.2)
# impute                 1.72.3    2023-01-19 [2] Bioconductor
# IRanges              * 2.32.0    2022-11-01 [2] Bioconductor
# irlba                  2.3.5.1   2022-10-03 [2] CRAN (R 4.2.1)
# iterators              1.0.14    2022-02-05 [2] CRAN (R 4.2.1)
# jaffelab             * 0.99.32   2023-02-09 [1] Github (LieberInstitute/jaffelab@7b7afe3)
# jsonlite               1.8.5     2023-06-05 [1] CRAN (R 4.2.3)
# KEGGREST               1.38.0    2022-11-01 [2] Bioconductor
# knitr                  1.42      2023-01-25 [2] CRAN (R 4.2.2)
# lattice                0.20-45   2021-09-22 [3] CRAN (R 4.2.3)
# lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.2.1)
# limma                * 3.54.2    2023-02-28 [2] Bioconductor
# locfit                 1.5-9.7   2023-01-02 [2] CRAN (R 4.2.2)
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.2.1)
# MASS                   7.3-58.2  2023-01-23 [3] CRAN (R 4.2.3)
# Matrix                 1.5-4     2023-04-04 [2] CRAN (R 4.2.3)
# MatrixGenerics       * 1.10.0    2022-11-01 [2] Bioconductor
# matrixStats          * 1.0.0     2023-06-02 [1] CRAN (R 4.2.3)
# memoise                2.0.1     2021-11-26 [2] CRAN (R 4.2.1)
# munsell                0.5.0     2018-06-12 [2] CRAN (R 4.2.1)
# nlme                   3.1-162   2023-01-31 [2] CRAN (R 4.2.2)
# nnet                   7.3-18    2022-09-28 [3] CRAN (R 4.2.3)
# pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.2.3)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.2.1)
# plyr                   1.8.8     2022-11-11 [2] CRAN (R 4.2.2)
# png                    0.1-8     2022-11-29 [2] CRAN (R 4.2.2)
# preprocessCore         1.60.2    2023-01-19 [2] Bioconductor
# prettyunits            1.1.1     2020-01-24 [2] CRAN (R 4.2.1)
# progress               1.2.2     2019-05-16 [2] CRAN (R 4.2.1)
# purrr                  1.0.1     2023-01-10 [2] CRAN (R 4.2.2)
# qvalue                 2.30.0    2022-11-01 [2] Bioconductor
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.2.1)
# rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 4.2.2)
# rappdirs               0.3.3     2021-01-31 [2] CRAN (R 4.2.1)
# RColorBrewer           1.1-3     2022-04-03 [2] CRAN (R 4.2.1)
# Rcpp                   1.0.10    2023-01-22 [2] CRAN (R 4.2.2)
# RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.2.3)
# readr                  2.1.4     2023-02-10 [2] CRAN (R 4.2.2)
# recount              * 1.24.1    2023-02-21 [2] Bioconductor
# rentrez                1.2.3     2020-11-10 [2] CRAN (R 4.2.1)
# reshape2               1.4.4     2020-04-09 [2] CRAN (R 4.2.1)
# restfulr               0.0.15    2022-06-16 [2] CRAN (R 4.2.1)
# rjson                  0.2.21    2022-01-09 [2] CRAN (R 4.2.1)
# rlang                  1.1.1     2023-04-28 [1] CRAN (R 4.2.3)
# rmarkdown              2.21      2023-03-26 [2] CRAN (R 4.2.3)
# rngtools               1.5.2     2021-09-20 [2] CRAN (R 4.2.1)
# rpart                  4.1.19    2022-10-21 [3] CRAN (R 4.2.3)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.2.1)
# Rsamtools              2.14.0    2022-11-01 [2] Bioconductor
# RSQLite                2.3.1     2023-04-03 [2] CRAN (R 4.2.3)
# rstudioapi             0.14      2022-08-22 [2] CRAN (R 4.2.1)
# rsvd                   1.0.5     2021-04-16 [2] CRAN (R 4.2.1)
# rtracklayer            1.58.0    2022-11-01 [2] Bioconductor
# S4Vectors            * 0.36.2    2023-02-26 [2] Bioconductor
# ScaledMatrix           1.6.0     2022-11-01 [2] Bioconductor
# scales                 1.2.1     2022-08-20 [2] CRAN (R 4.2.1)
# scater               * 1.26.1    2022-11-13 [2] Bioconductor
# scuttle              * 1.8.4     2023-01-19 [2] Bioconductor
# segmented              1.6-4     2023-04-13 [1] CRAN (R 4.2.3)
# sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.2.1)
# SingleCellExperiment * 1.20.1    2023-03-17 [2] Bioconductor
# sparseMatrixStats      1.10.0    2022-11-01 [2] Bioconductor
# stringi                1.7.12    2023-01-11 [2] CRAN (R 4.2.2)
# stringr                1.5.0     2022-12-02 [2] CRAN (R 4.2.2)
# SummarizedExperiment * 1.28.0    2022-11-01 [2] Bioconductor
# survival               3.5-3     2023-02-12 [3] CRAN (R 4.2.3)
# tibble                 3.2.1     2023-03-20 [2] CRAN (R 4.2.3)
# tidyr                  1.3.0     2023-01-24 [2] CRAN (R 4.2.2)
# tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.2.1)
# tzdb                   0.3.0     2022-03-28 [2] CRAN (R 4.2.1)
# utf8                   1.2.3     2023-01-31 [2] CRAN (R 4.2.2)
# VariantAnnotation      1.44.1    2023-02-15 [2] Bioconductor
# vctrs                  0.6.2     2023-04-19 [1] CRAN (R 4.2.3)
# vipor                  0.4.5     2017-03-22 [2] CRAN (R 4.2.1)
# viridis                0.6.2     2021-10-13 [2] CRAN (R 4.2.1)
# viridisLite            0.4.2     2023-05-02 [1] CRAN (R 4.2.3)
# WGCNA                * 1.72-1    2023-01-18 [1] CRAN (R 4.2.2)
# withr                  2.5.0     2022-03-03 [2] CRAN (R 4.2.1)
# xfun                   0.39      2023-04-20 [1] CRAN (R 4.2.3)
# XML                    3.99-0.14 2023-03-19 [2] CRAN (R 4.2.3)
# xml2                   1.3.3     2021-11-30 [2] CRAN (R 4.2.1)
# XVector                0.38.0    2022-11-01 [2] Bioconductor
# yaml                   2.3.7     2023-01-23 [2] CRAN (R 4.2.2)
# zlibbioc               1.44.0    2022-11-01 [2] Bioconductor
# 
# [1] /users/bsimbiat/R/4.2.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/library
# 
# ─────────────────────────────────────────────────────────────────────────────────────
