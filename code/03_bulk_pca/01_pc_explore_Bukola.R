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


# Before Brain Swaps ###########################################################
# Loading relevant rse objects #################################################
# gene
load(here("processed-data", "02_bulk_qc", "count_data_bukola", 
          "rse_gene_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata"))
rse_gene = rse_gene_filt
rm(rse_gene_filt)

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
pc_rse_gene[[1]]$"AgeDeath" = as.numeric(pc_rse_gene[[1]]$"AgeDeath")

## Creating Plotting PC base function ##########################################

pc_to_pc <- function (pcx, pcy, pc_df, colorby) {
  
  # unlisting returned object from pca_creator function
  pc_data = pc_df[[1]]
  pc_variables = pc_df[[2]]
  
  # Use pos to ensure jitter and text_repel share coordinates (prevents mislabeling).  
  pos <- position_jitter(seed = 2)
  
  # grabbing titles 
  x_titler = pc_variables[grepl(paste0(pcx, ":"), pc_variables)]
  y_titler = pc_variables[grepl(paste0(pcy, ":"), pc_variables)]
  
  # plotting with Brain Numbers
  # for discrete colorby (discrete phenotypes)

# if (is.factor(pc_data[, colorby]) == TRUE){
  plot = ggplot(pc_data, aes_(x = pc_data[, pcx], y = pc_data[, pcy])) + 
    geom_jitter(aes_(color = pc_data[,colorby]), position = pos) +
    geom_text_repel(aes(label = pc_data[,"BrNum"], color = 
                          pc_data[,colorby]), position = pos) +
    labs(x = x_titler, y = y_titler) +
    theme_bw(base_size = 10) + 
    theme(legend.position= "bottom", plot.margin=unit (c (1.5,2,1,2), 'cm'), 
          axis.text.x = element_text(vjust = 0.7), text = element_text(size=15),
          axis.title = element_text(size=15)) +
    guides(color = guide_legend(title =  rename_vars[rename_vars$orig_var_name == 
           colorby,]$var_plot_title))
 #  } else {
 #    
 #   plot = ggplot(pc_data, aes_(x = pc_data[, pcx], y = pc_data[, pcy]) +
 #     scale_colour_gradient2(colours = pc_data[,colorby])) +
 #     geom_jitter(aes_(color = pc_data[,colorby]), position = pos) +
 #     geom_text_repel(aes(label = pc_data[,"BrNum"]), position = pos) +
 #     labs(x = x_titler, y = y_titler) +
 #     theme_bw(base_size = 10) + 
 #     theme(legend.position= "bottom", plot.margin=unit (c (1.5,2,1,2), 'cm'), 
 #           axis.text.x = element_text(vjust = 0.7), text = element_text(size=15),
 #          axis.title = element_text(size=15)) +
 #     guides(color = guide_legend(title =  rename_vars[rename_vars$orig_var_name == 
 #                                                        colorby,]$var_plot_title))
 # }
  
  return(plot)
}


## Plot and save

# PC1 vs PC2
pdf(file = here("plots", "03_bulk_pca", "pc_plots_bukola", "PC1_vs_PC2_rse_gene_n69.pdf"))
  lapply(phenoCols, FUN = pc_to_pc, pcx = "PC1", pcy = "PC2", pc_df = pc_rse_gene)
dev.off()

# PC2 vs PC3
pdf(file = here("plots", "03_bulk_pca", "pc_plots_bukola", "PC2_vs_PC3_rse_gene_n69.pdf"))
lapply(phenoCols, FUN = pc_to_pc, pcx = "PC2", pcy = "PC3", pc_df = pc_rse_gene)
dev.off()

# PC3 vs PC4
pdf(file = here("plots", "03_bulk_pca", "pc_plots_bukola", "PC3_vs_PC4_rse_gene_n69.pdf"))
lapply(phenoCols, FUN = pc_to_pc, pcx = "PC3", pcy = "PC4", pc_df = pc_rse_gene)
dev.off()

# PC4 vs PC5
pdf(file = here("plots", "03_bulk_pca", "pc_plots_bukola", "PC4_vs_PC5_rse_gene_n69.pdf"))
lapply(phenoCols, FUN = pc_to_pc, pcx = "PC4", pcy = "PC5", pc_df = pc_rse_gene)
dev.off()


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
