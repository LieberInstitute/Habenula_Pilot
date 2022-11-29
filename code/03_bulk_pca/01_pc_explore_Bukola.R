# November 29, 2022
# 01_pc_explore_Bukola.R - Creating  PCA plots for bulk RNA-seq data pre and post
# brain sample drops.
# qrsh -l mem_free=20G,h_vmem=20G

library(SummarizedExperiment)
library(recount)
library(jaffelab)
library(RColorBrewer)
library(here)

# Before Brain Swaps ###########################################################
# Loading relevant rse objects #################################################
# gene
load(here("processed-data", "02_bulk_qc", "count_data_bukola", 
          "rse_gene_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata"))
rse_gene = rse_gene_filt

# exon: NEED TO REMAKE
# load(here("processed-data", "02_bulk_qc", "count_data_bukola", 
#           "rse_exon_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata")) 
# rse_exon = rse_exon_filt
 
# jx
load(here("processed-data", "02_bulk_qc", "count_data_bukola", 
          "rse_jx_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata")) 
rse_jx = rse_jx_filt
 
# tx
load(here("processed-data", "02_bulk_qc", "count_data_bukola", 
          "rse_tx_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata"))
rse_tx = rse_tx_filt

## Creating PCA base function:
PCA <- function(rse_type){
  
  pca <- prcomp(t(assays(rse_type)$logcounts))
  
  ## % of the variance explained by each PC
  pca_vars <- getPcaVars(pca)
  pca_vars_labs<- paste0("PC", seq(along = pca_vars), ": ",
    pca_vars, "% Var Expl")
  
  ## Joining PC and sample info
  pca_data<-cbind(pca$x, colData(rse_type))
  
  ## Add samples' phenotypes
  pca_data<-as.data.frame(pca_data)
  
  # Returning PCA data as well as their plottable labels.
  return(list(pca_data, pca_vars_labs))
}

## Applying PCA function per rse type
pc_rse_gene = PCA(rse_gene)
pc_rse_gene = PCA(rse_gene)
pc_rse_gene = PCA(rse_gene)
pc_rse_gene = PCA(rse_gene)

## Plotting with geom_jitter

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
