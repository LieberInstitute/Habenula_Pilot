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

## Coloring by
colorbylist = c("PrimaryDx", "Flowcell", "AgeDeath", "RIN", "percentGC_R1",
                "percentGC_R2", "ERCCsumLogErr", "overallMapRate", "concordMapRate",
                "totalMapped", "mitoMapped", "mitoRate", "totalAssignedGene", 
                "rRNA_rate", "detected", "logNumReads", "logNumMapped",
                "logNumUnmapped")

## FUNCTION 1: CREATES PCA VALUES FOR RSE OBJECT ###############################
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

## FUNCTION 2: CREATES PCA PLOTS ###############################################
pc_to_pc <- function (pcx, pcy, pc_df, colorbylist, dataType, numdrop = NA, 
                      sampsdropped = NA) {
  # pcx and pcy are the pcas of interest
  # pc_df is the output object from the PCA creator function
  # colorbylist is the list of metrics to color plots buy (changes number of 
  # plots per saved file)
  # dataType options are: "exon" "tx" "jx" "gene"
  # sampsdropped is one character string of brains that we dropped 
  # (numdrop and sampsdropped are either both NA OR both not)
  
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

     if (i == "PrimaryDx"){
      
       forhighlight = pc_data[pc_data$BrNum %in% c("Br1676","Br5459", "Br6323"),]
       nothighlight = pc_data[! pc_data$BrNum %in% c("Br1676","Br5459", "Br6323"),]
         
       plot = ggplot(pc_data, aes_string(x = pcx, y = pcy)) +
        geom_jitter(aes_string(colour = i), data = nothighlight, alpha = 0.8, position = pos, size = 5) +
         geom_jitter(aes_string(x = pcx, y = pcy), colour = "red", 
                     data = forhighlight,alpha = 0.8, position = pos, size = 5) +
        geom_text_repel(aes_string(label = "BrNum"), position = pos, 
                        color = "lightgrey", max.overlaps = 5) +
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
          geom_text_repel(aes_string(label = "BrNum"), position = pos, 
                          color = "lightgrey") +
          labs (x = x_titler, y = y_titler, 
                color = rename_vars[rename_vars$orig_var_name ==  
                                      i,]$var_plot_title) +
          theme_bw(base_size = 10) +
          theme(legend.position= "bottom", legend.key.width = unit(1.5, 'cm'),
                plot.margin=unit (c (1.5,2,1,2), 'cm'),
                axis.text.x = element_text(vjust = 0.7), 
                text = element_text(size=15),
                axis.title = element_text(size=15))  
       
      }
    
    plot_list[[c]] = plot 
    c = c + 1
}
  # Saving plots is easier in the function
  firstnamer = paste(pcx, "vs", pcy, "rse", dataType, sep = "_")
  secondnamer = paste0("_n", as.character(dim(pc_data)[1]), ".pdf") 
  nameFILE = paste0(firstnamer, secondnamer)
  
  if (is.na(numdrop) == TRUE) {
    type = "before_drop"
  
    pdf(file = here("plots", "03_bulk_pca", "pc_plots_bukola", type, nameFILE))
      print(plot_list)
    dev.off()
    
  } else if(is.na(numdrop) == FALSE && is.na(sampsdropped) == TRUE) { 
    
    print("Error. Which samples are we dropping?")
    
  } else if(is.na(numdrop) == FALSE && is.na(sampsdropped) == FALSE){
    
    type = paste0("dropped", numdrop, "brain")
    orgbybrain = paste0("pc_plots_wout_", sampsdropped)
    
    dir.create(here("plots", "03_bulk_pca", "pc_plots_bukola", type, orgbybrain))
    
    pdf(file = here("plots", "03_bulk_pca", "pc_plots_bukola", type, orgbybrain, nameFILE))
      print(plot_list)
    dev.off()
  }
  
}

## FUNCTION 3: TO FROP SAMPLES #################################################
# Creating dropping function
dropSample <- function(rse_obj, sampler) {
  # rse_obj is the rse object we would like to drop from (specific to rse_gene)
  # sample should be character values in list form :)
  sampler = as.list(sampler)
  
  pd = colData(rse_obj)
  finder = pd[pd$BrNum %in% sampler, ]$RNum
  
  pd = pd[! pd$RNum %in% finder,]
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

### FUNCTION 4: MASTER FUNCTION TO RUN ALL RELEVANT FUNCTIONS ABOVE ############
# MASTER: Most upstream function to fun all relevant functions on rse obj of interst
pca_print <- function(rse_obj, colorbylist, dataType, numdrop, sampsdropped){
  # rse_obj is output from dropSample or regular rse obj
  # colorbylist is metrics to color plots by (list)
  # dataType has to be "gene", "exon", "tx", or "jx"
  # numdrop 
  
  # Calc PC:
  pc_rse = pca_creator(rse_obj)
  
  # Updating variable types 
  pc_rse[[1]]$"PrimaryDx" = as.factor(pc_rse[[1]]$"PrimaryDx")
  pc_rse[[1]]$"Flowcell" = as.factor(pc_rse[[1]]$"Flowcell")
  pc_rse[[1]]$"percentGC_R1" = as.factor(pc_rse[[1]]$"percentGC_R1")
  pc_rse[[1]]$"percentGC_R2" = as.factor(pc_rse[[1]]$"percentGC_R2")
  
  
  ## PC1 vs PC2
  pc_to_pc("PC1", "PC2", pc_df = pc_rse, colorbylist, dataType, numdrop, sampsdropped)
  ## PC2 vs PC3
  pc_to_pc("PC2", "PC3", pc_df = pc_rse, colorbylist, dataType, numdrop, sampsdropped)
  ## PC3 vs PC4
  pc_to_pc("PC3", "PC4", pc_df = pc_rse, colorbylist, dataType, numdrop, sampsdropped)
  ## PC4 vs PC5
  pc_to_pc("PC4", "PC5", pc_df = pc_rse, colorbylist, dataType, numdrop, sampsdropped)

}

### RUNNING CODE ON DATA #######################################################
## BEFORE DROPS ################################################################ 
## Gene:
pc_rse_gene = pca_creator(rse_gene) # gene only for now 

## Correcting rse objects 
colData(rse_gene)$PrimaryDx <- recode_factor(colData(rse_gene)$PrimaryDx, 
                                             Schizo = "SCZD")
# Adding log 10 of nums 
colData(rse_gene)$logNumReads <- log10(colData(rse_gene)$numReads)
colData(rse_gene)$logNumMapped <- log10(colData(rse_gene)$numMapped)
colData(rse_gene)$logNumUnmapped <- log10(colData(rse_gene)$numUnmapped)

# Running all codes with master:
pca_print(rse_gene, colorbylist, "gene", numdrop = NA, sampsdropped = NA)
  

# DROPPING 1 SAMPLE ############################################################
rse_drop1676 = dropSample(rse_gene, c("Br1676"))
rse_drop5459 = dropSample(rse_gene, c("Br5459"))
rse_drop6323 = dropSample(rse_gene, c("Br6323"))

########## PCA after dropping Br1676 ###########################################
colData(rse_drop1676)$PrimaryDx <- recode_factor(colData(rse_drop1676)$PrimaryDx, 
                                             Schizo = "SCZD")
# Adding log 10 of nums 
colData(rse_drop1676)$logNumReads <- log10(colData(rse_drop1676)$numReads)
colData(rse_drop1676)$logNumMapped <- log10(colData(rse_drop1676)$numMapped)
colData(rse_drop1676)$logNumUnmapped <- log10(colData(rse_drop1676)$numUnmapped)

# Running all codes with master:
pca_print(rse_drop1676, colorbylist, "gene", numdrop = 1, sampsdropped = "Br1676")

########## PCA after dropping Br5459 ###########################################
colData(rse_drop5459)$PrimaryDx <- recode_factor(colData(rse_drop5459)$PrimaryDx, 
                                                 Schizo = "SCZD")
# Adding log 10 of nums 
colData(rse_drop5459)$logNumReads <- log10(colData(rse_drop5459)$numReads)
colData(rse_drop5459)$logNumMapped <- log10(colData(rse_drop5459)$numMapped)
colData(rse_drop5459)$logNumUnmapped <- log10(colData(rse_drop5459)$numUnmapped)

# Running all codes with master:
pca_print(rse_drop5459, colorbylist, "gene", numdrop = 1, sampsdropped = "Br5459")

########## PCA after dropping Br6323 ###########################################
colData(rse_drop6323)$PrimaryDx <- recode_factor(colData(rse_drop6323)$PrimaryDx, 
                                                 Schizo = "SCZD")
# Adding log 10 of nums 
colData(rse_drop6323)$logNumReads <- log10(colData(rse_drop6323)$numReads)
colData(rse_drop6323)$logNumMapped <- log10(colData(rse_drop6323)$numMapped)
colData(rse_drop6323)$logNumUnmapped <- log10(colData(rse_drop6323)$numUnmapped)

# Running all codes with master:
pca_print(rse_drop6323, colorbylist, "gene", numdrop = 1, sampsdropped = "Br6323")

# DROPPING 2 SAMPLES ###########################################################
rse_drop1676and5459 = dropSample(rse_gene, c("Br1676", "Br5459"))
rse_drop5459and6323 = dropSample(rse_gene, c("Br5459", "Br6323"))

########## PCA after dropping Br1676 and Br5459 ################################
colData(rse_drop1676and5459)$PrimaryDx <- 
  recode_factor(colData(rse_drop1676and5459)$PrimaryDx, Schizo = "SCZD")

# Adding log 10 of nums 
colData(rse_drop1676and5459)$logNumReads <- log10(colData(rse_drop1676and5459)$numReads)
colData(rse_drop1676and5459)$logNumMapped <- log10(colData(rse_drop1676and5459)$numMapped)
colData(rse_drop1676and5459)$logNumUnmapped <- log10(colData(rse_drop1676and5459)$numUnmapped)

# Running all codes with master:
pca_print(rse_drop1676and5459, colorbylist, "gene", numdrop = 2, 
          sampsdropped = c("Br1676_Br5459"))

########## PCA after dropping Br5459 and Br6323 ################################
colData(rse_drop5459and6323)$PrimaryDx <- 
  recode_factor(colData(rse_drop5459and6323)$PrimaryDx, Schizo = "SCZD")

# Adding log 10 of nums 
colData(rse_drop5459and6323)$logNumReads <- log10(colData(rse_drop5459and6323)$numReads)
colData(rse_drop5459and6323)$logNumMapped <- log10(colData(rse_drop5459and6323)$numMapped)
colData(rse_drop5459and6323)$logNumUnmapped <- log10(colData(rse_drop5459and6323)$numUnmapped)

# Running all codes with master:
pca_print(rse_drop5459and6323, colorbylist, "gene", numdrop = 2, 
          sampsdropped = c("Br5459_Br6323"))

# DROPPING ALL 3 SAMPLES #######################################################
rse_drop_all3 = dropSample(rse_gene, c("Br5459", "Br6323", "Br1676"))

########## PCA after dropping all 3 ############################################
colData(rse_drop_all3)$PrimaryDx <- 
  recode_factor(colData(rse_drop_all3)$PrimaryDx, Schizo = "SCZD")

# Adding log 10 of nums 
colData(rse_drop_all3)$logNumReads <- log10(colData(rse_drop_all3)$numReads)
colData(rse_drop_all3)$logNumMapped <- log10(colData(rse_drop_all3)$numMapped)
colData(rse_drop_all3)$logNumUnmapped <- log10(colData(rse_drop_all3)$numUnmapped)

# Running all codes with master:
pca_print(rse_drop_all3, colorbylist, "gene", numdrop = 2, 
          sampsdropped = c("all3"))


### printing tester
pdf(file = here("plots", "03_bulk_pca", "pc_plots_bukola", "before_drop", "tester.pdf"))
  plot 
dev.off()

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
