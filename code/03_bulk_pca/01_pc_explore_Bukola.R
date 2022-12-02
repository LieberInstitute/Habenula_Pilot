# November 29, 2022
# 01_pc_explore_Bukola.R - Creating  PCA plots for bulk RNA-seq data pre and post
# brain sample drops.
# qrsh -l mem_free=20G,h_vmem=20G

library(SummarizedExperiment)
library(recount)
library(jaffelab)
library(RColorBrewer)
library(here)
library(ggrepel)
library(viridis)
library(dplyr)
library(VariantAnnotation)
library(WGCNA) 
library(biomartr) 
library(edgeR)
library(sessioninfo)


# Before Brain Swaps ###########################################################
# Loading relevant rse objects #################################################
# gene
load(here("processed-data", "02_bulk_qc", "count_data_bukola", 
          "rse_gene_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata"))

# exon
load(here("processed-data", "02_bulk_qc", "count_data_bukola", 
         "rse_exon_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata")) 

# jx
load(here("processed-data", "02_bulk_qc", "count_data_bukola", 
          "rse_jx_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata")) 

# tx
load(here("processed-data", "02_bulk_qc", "count_data_bukola", 
          "rse_tx_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata"))

## Correcting rse objects accordingly ##########################################
# Changing schizo to SCZD
colData(rse_gene)$PrimaryDx <- recode_factor(colData(rse_gene)$PrimaryDx, Schizo = "SCZD")
# Adding log 10 of nums 
colData(rse_gene)$logNumReads <- log10(colData(rse_gene)$numReads)
colData(rse_gene)$logNumMapped <- log10(colData(rse_gene)$numMapped)
colData(rse_gene)$logNumUnmapped <- log10(colData(rse_gene)$numUnmapped)


## Including relevant stable variables #########################################
# Variable phenotypes in our data 

phenoCols = as.vector(c("PrimaryDx", "Flowcell", "AgeDeath"))

# Relevant QC metrics  
QCmetCols = c("RIN", "percentGC_R1", "percentGC_R2", "ERCCsumLogErr",
              "overallMapRate", "concordMapRate", "totalMapped", 
              "mitoMapped", "mitoRate", "totalAssignedGene", "rRNA_rate", 
              "detected", "logNumReads", "logNumMapped", "logNumUnmapped")

# rename_vars ####
# Creating df for plot text to rename variables:
orig_var_name <- c("RNum", "RIN", "BrNum", "AgeDeath", "Sex", "PrimaryDx", 
                   "percentGC_R1", "percentGC_R2", "ERCCsumLogErr", 
                   "numReads", "numMapped", "numUnmapped", "overallMapRate", 
                   "concordMapRate", "totalMapped", "mitoMapped", "mitoRate", 
                   "totalAssignedGene", "rRNA_rate", "Flowcell", 
                   "sum", "detected", "subsets_Mito_sum", "subsets_Mito_detected",
                   "subsets_Mito_percent", "subsets_Ribo_sum", 
                   "subsets_Ribo_detected", "subsets_Ribo_percent", "logNumReads",
                   "logNumMapped", "logNumUnmapped", "AgeInterval")

var_plot_title <- c("RNum", "RIN", "Brain Number", "Age of Death", "Sex",
                    "Primary Dx", "Percent GC R1", "Percent GC R2", 
                    "ERCC RSS", "Num of Reads", "Num Mapped", "Num Unmapped", 
                    "Overall Map Rate", "Concordant Map Rate", "Total Mapped", 
                    "chrM Mapped", "chrM Map Rate", "Total Assigned Genes", 
                    "Gene rRNA Rate", "Flowcell", "Sum", "Detected", 
                    "Sum of Mito Subsets", "Detected Mito Subsets", 
                    "Percent of Mito Subsets", "Sum of Ribo Subsets", 
                    "Detected Ribo Subsets", "Percent of Ribo Subsets",
                    "Num Reads (log 10)", "Num Mapped (log 10)", 
                    "Num Unmapped (log 10)", "Age Intervals")

rename_vars <- data.frame(orig_var_name, var_plot_title)

## Creating PCA base function ##################################################
pca_creator <- function(rse_type){
  
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

## Applying PCA function per rse type ##########################################
pc_rse_gene = pca_creator(rse_gene) # gene only for now 

## Updating the variable type for phenotype (discrete vs continuous) ###########
pc_rse_gene[[1]]$"PrimaryDx" = as.factor(pc_rse_gene[[1]]$"PrimaryDx")
pc_rse_gene[[1]]$"Flowcell" = as.factor(pc_rse_gene[[1]]$"Flowcell")
pc_rse_gene[[1]]$"percentGC_R1" = as.factor(pc_rse_gene[[1]]$"percentGC_R1")
pc_rse_gene[[1]]$"percentGC_R2" = as.factor(pc_rse_gene[[1]]$"percentGC_R2")


## Creating Plotting PC base function ##########################################
pc_to_pc <- function (pcx, pcy, pc_df, colorbylist, dataType, numdrop = NA) {
  # pcx and pcy are the pcas of interest
  # pc_df is the output object from the PCA creator function
  # colorbylist is the list of metrics to color plots buy (changes number of 
  # plots per saved file)
  # dataType options are: "exon" "tx" "jx" "gene"
  # folder = "before_drop"
  
  # unlisting returned object from pca_creator function
  pc_data = pc_df[[1]]
  pc_variables = pc_df[[2]] 
  
  # Use pos to ensure jitter and text_repel share coordinates (prevents mislabeling).  
  pos <- position_jitter(seed = 2)
  
  # grabbing titles 
  x_titler = pc_variables[grepl(paste0(pcx, ":"), pc_variables)]
  y_titler = pc_variables[grepl(paste0(pcy, ":"), pc_variables)]
  
  # prepping for saving plots into list
  plot_list = list()
  c = 1
 
  # plotting by coloring scheme
  for(i in colorbylist){

     if (is.factor(pc_data[, i]) == TRUE){
      
       plot = ggplot(pc_data, aes_string(x = pcx, y = pcy)) +
        scale_color_viridis(option = "plasma", discrete = TRUE) +
        geom_jitter(aes_string(colour = i), position = pos, size = 7) +
        geom_text_repel(aes_string(label = "BrNum", color = i), position = pos) +
        labs(x = x_titler, y = y_titler,
             color = rename_vars[rename_vars$orig_var_name ==  i,]$var_plot_title) +
        theme_bw(base_size = 10) + 
        theme(legend.position= "bottom", plot.margin=unit (c (1.5,2,1,2), 'cm'), 
              axis.text.x = element_text(vjust = 0.7), text = element_text(size=15),
              axis.title = element_text(size=15))
       
      } else {

        plot = ggplot(pc_data, aes_string(x = pcx, y = pcy)) +
          scale_colour_viridis() +
          geom_jitter(aes_string(colour = i), position = pos, size = 7) +
          geom_text_repel(aes_string(label = "BrNum"), position = pos, color = "lightgrey") +
          labs (x = x_titler, y = y_titler, 
                color = rename_vars[rename_vars$orig_var_name ==  i,]$var_plot_title) +
          theme_bw(base_size = 10) +
          theme(legend.position= "bottom", legend.key.width = unit(1.5, 'cm'),
                plot.margin=unit (c (1.5,2,1,2), 'cm'),
                axis.text.x = element_text(vjust = 0.7), text = element_text(size=15),
                axis.title = element_text(size=15))  
       
      }
    
    plot_list[[c]] = plot 
    c = c + 1
}
  
  if (is.na(numdrop) == TRUE) {
    type = "before_drop"
  } else { 
    type = paste0("dropped", numdrop, "brain")
  }
  
  # Saving plots is easier in the function
  firstnamer = paste(pcx, "vs", pcy, "rse", dataType, sep = "_")
  secondnamer = paste0("_n", as.character(dim(pc_data)[1]), ".pdf") 
  nameFILE = paste0(firstnamer, secondnamer)
  
  pdf(file = here("plots", "03_bulk_pca", "pc_plots_bukola", type, nameFILE))
    print(plot_list)
  dev.off()
}

## Coloring by

colorbylist = c("PrimaryDx", "Flowcell", "AgeDeath", "RIN", "percentGC_R1",
                "percentGC_R2", "ERCCsumLogErr", "overallMapRate", "concordMapRate",
                "totalMapped", "mitoMapped", "mitoRate", "totalAssignedGene", 
                "rRNA_rate", "detected", "logNumReads", "logNumMapped",
                "logNumUnmapped")

## BEFORE DROP #################################################################
## rse_gene ####################################################################
## PC1 vs PC2
pc_to_pc("PC1", "PC2", pc_df = pc_rse_gene, colorbylist, dataType = "gene")
## PC2 vs PC3
pc_to_pc("PC2", "PC3", pc_df = pc_rse_gene, colorbylist, dataType = "gene")
## PC3 vs PC4
pc_to_pc("PC3", "PC4", pc_df = pc_rse_gene, colorbylist, dataType = "gene")
## PC4 vs PC5
pc_to_pc("PC4", "PC5", pc_df = pc_rse_gene, colorbylist, dataType = "gene")


## DROPPING SAMPLES ############################################################
## Fixing RNums onto colnames of rse
colnames(rse_gene) <- rse_gene$RNum

# Creating dropping function
dropSample <- function(rse_obj, sample) {
  # rse_obj is the rse object we would like to drop from (specific to rse_gene)
  # sample should be character values in list form :)
  pd = colData(rse_obj)
  finder = pd[pd$BrNum == sample, ]$RNum
  
  pd = pd[pd$RNum != finder,]
  assays(rse_obj) <- assays(rse_obj)[1]
  
  rse_obj <- rse_obj[, pd$RNum]
  
  testInteg = sum(rowSums(is.na(assay(rse_obj)) | assay(rse_obj) == "") > 0) 
  
  if(testInteg != 0){
    return("Error. NAs introduced after drop.")
  }
  
  determine = (length(which(assay(rse_obj) == 0)) * 100) / 
    (nrow(rse_obj)*ncol(rse_obj))

  if (determine > 0.75){
    assays(rse_obj, withDimnames=FALSE)$logcounts<- 
      edgeR::cpm(calcNormFactors(rse_obj, method = "TMMwsp"), 
                 log = TRUE, prior.count = 0.5)
  } else if (determine <= 0.75){
    assays(rse_obj, withDimnames=FALSE)$logcounts = 
      edgeR::cpm(calcNormFactors(rse_obj, method = "TMM"), 
                 log = TRUE, prior.count = 0.5)
  }

  return(rse_obj)
}

## Testing drop function

rse_gene_drop1a = dropSample(rse_gene, c("Br1676"))


# Drop 1: 
rse_gene_drop1 = rse_gene
pd = colData(rse_gene_drop1) 

# finding RNum of brain sample we want to drop
drop1676 = pd[pd$BrNum == "Br1676", ]$RNum

# dropping sample and re-saving gene counts
pd = pd[pd$RNum != drop1676,]
assays(rse_gene_drop1) <- assays(rse_gene_drop1)[1]
rse_gene_drop1 <- rse_gene_drop1[, pd$RNum]

# verify data integrity 
which(rowSums(is.na(assay(rse_gene_drop1)) | assay(rse_gene_drop1) == "") > 0) 
# none

(length(which(assay(rse_gene_drop1) == 0)) * 100) / (nrow(rse_gene_drop1)*ncol(rse_gene_drop1))
# 0.75 which is high proportion of zeros so TMMwsp method instead of TMM for logcounts

# calc logcounts
assays(rse_gene_drop1, withDimnames=FALSE)$logcounts<- 
  edgeR::cpm(calcNormFactors(rse_gene_drop1, method = "TMMwsp"), 
             log = TRUE, prior.count = 0.5)

# Calc PC:
pc_rse_gene_drop1 = pca_creator(rse_gene_drop1)

# Updating variable types 
pc_rse_gene_drop1[[1]]$"PrimaryDx" = as.factor(pc_rse_gene_drop1[[1]]$"PrimaryDx")
pc_rse_gene_drop1[[1]]$"Flowcell" = as.factor(pc_rse_gene_drop1[[1]]$"Flowcell")
pc_rse_gene_drop1[[1]]$"percentGC_R1" = as.factor(pc_rse_gene_drop1[[1]]$"percentGC_R1")
pc_rse_gene_drop1[[1]]$"percentGC_R2" = as.factor(pc_rse_gene_drop1[[1]]$"percentGC_R2")


# Plotting after drop of Br1676
## PC1 vs PC2
pc_to_pc("PC1", "PC2", pc_df = pc_rse_gene_drop1, colorbylist, dataType = "gene", numdrop = 1)
## PC2 vs PC3
pc_to_pc("PC2", "PC3", pc_df = pc_rse_gene_drop1, colorbylist, dataType = "gene", numdrop = 1)
## PC3 vs PC4
pc_to_pc("PC3", "PC4", pc_df = pc_rse_gene_drop1, colorbylist, dataType = "gene", numdrop = 1)
## PC4 vs PC5
pc_to_pc("PC4", "PC5", pc_df = pc_rse_gene_drop1, colorbylist, dataType = "gene", numdrop = 1)



## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
