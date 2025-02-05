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
library(wesanderson)


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
    
    forhighlight = pc_data[pc_data$BrNum %in% c("Br1676","Br5459", "Br6323"),]
    nothighlight = pc_data[! pc_data$BrNum %in% c("Br1676","Br5459", "Br6323"),]
    
     if (i == "PrimaryDx"){
       plot = ggplot(pc_data, aes_string(x = pcx, y = pcy)) +
        scale_color_manual(values=wes_palette(n=2, name="GrandBudapest1"))+
        geom_jitter(aes_string(colour = i), data = nothighlight, alpha = 0.8, position = pos, size = 6) +
         geom_jitter(aes_string(x = pcx, y = pcy), colour = "darkblue", 
                     data = forhighlight, alpha = 0.8, position = pos, size = 6) +
        geom_text_repel(aes_string(label = "BrNum"), position = pos, 
                        color = "lightgrey", max.overlaps = 5) +
        labs(x = x_titler, y = y_titler,
             color = rename_vars[rename_vars$orig_var_name ==  i,]$var_plot_title,
             caption = "Br1676, Br5459, and Br6323") +
        theme_bw(base_size = 10) + 
        theme(legend.position= "bottom", plot.margin=unit (c (1.5,2,1,2), 'cm'), 
              axis.text.x = element_text(vjust = 0.7), text = element_text(size=15),
              axis.title = element_text(size=15),
              plot.caption = element_text(color = "darkblue", face = "italic")) 
       
     } else if (i == "Flowcell"){
       
       plot = ggplot(pc_data, aes_string(x = pcx, y = pcy)) +
         scale_color_manual(values=wes_palette(n=2, name="IsleofDogs1"))+
         geom_jitter(aes_string(colour = i), data = nothighlight, alpha = 0.8, position = pos, size = 6) +
         geom_jitter(aes_string(x = pcx, y = pcy), colour = "darkblue", 
                     data = forhighlight, alpha = 0.8, position = pos, size = 6) +
         geom_text_repel(aes_string(label = "BrNum"), position = pos, 
                         color = "lightgrey", max.overlaps = 5) +
         labs(x = x_titler, y = y_titler,
              color = rename_vars[rename_vars$orig_var_name ==  i,]$var_plot_title,
              caption = "Br1676, Br5459, and Br6323") +
         theme_bw(base_size = 10) + 
         theme(legend.position= "bottom", plot.margin=unit (c (1.5,2,1,2), 'cm'), 
               axis.text.x = element_text(vjust = 0.7), text = element_text(size=15),
               axis.title = element_text(size=15),
               plot.caption = element_text(color = "darkblue", face = "italic")) 
       
      } else {

        plot = ggplot(pc_data, aes_string(x = pcx, y = pcy)) +
          scale_color_gradientn(colors = wes_palette("Darjeeling1", type = "continuous")) +
          geom_jitter(aes_string(colour = i), data = nothighlight, alpha = 0.7,
                      position = pos, size = 5) +
          geom_jitter(aes_string(x = pcx, y = pcy), colour = "darkblue", 
                      data = forhighlight, alpha = 0.8, position = pos, size = 6) +
          geom_text_repel(aes_string(label = "BrNum"), position = pos, 
                          color = "lightgrey", max.overlaps = 6) +
          labs(x = x_titler, y = y_titler,
               color = rename_vars[rename_vars$orig_var_name ==  i,]$var_plot_title,
               caption = "Br1676, Br5459, and Br6323") +
          theme_bw(base_size = 10) + 
          theme(legend.position= "bottom", plot.margin=unit (c (1.5,2,1,2), 'cm'), 
                axis.text.x = element_text(vjust = 0.7), text = element_text(size=15),
                axis.title = element_text(size=15),
                plot.caption = element_text(color = "darkblue", face = "italic")) 
       
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

# ─ Session info ─────────────────────────────────────────────────────────────────────────────────────────────────────────────
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
# ─ Packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# AnnotationDbi          1.60.2    2023-03-10 [2] Bioconductor
# backports              1.4.1     2021-12-13 [2] CRAN (R 4.2.1)
# base64enc              0.1-3     2015-07-28 [2] CRAN (R 4.2.1)
# Biobase              * 2.58.0    2022-11-01 [2] Bioconductor
# BiocFileCache          2.6.1     2023-02-17 [2] Bioconductor
# BiocGenerics         * 0.44.0    2022-11-01 [2] Bioconductor
# BiocIO                 1.8.0     2022-11-01 [2] Bioconductor
# BiocParallel           1.32.6    2023-03-17 [2] Bioconductor
# biomaRt                2.54.1    2023-03-20 [2] Bioconductor
# biomartr             * 1.0.3     2023-05-07 [1] CRAN (R 4.2.3)
# Biostrings           * 2.66.0    2022-11-01 [2] Bioconductor
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
# derfinder              1.32.0    2022-11-01 [2] Bioconductor
# derfinderHelper        1.32.0    2022-11-01 [2] Bioconductor
# digest                 0.6.31    2022-12-11 [2] CRAN (R 4.2.2)
# doParallel             1.0.17    2022-02-07 [2] CRAN (R 4.2.1)
# doRNG                  1.8.6     2023-01-16 [2] CRAN (R 4.2.2)
# downloader             0.4       2015-07-09 [2] CRAN (R 4.2.1)
# dplyr                * 1.1.1     2023-03-22 [2] CRAN (R 4.2.3)
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
# ggplot2              * 3.4.2     2023-04-03 [2] CRAN (R 4.2.3)
# ggrepel              * 0.9.3     2023-02-03 [2] CRAN (R 4.2.2)
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
# RColorBrewer         * 1.1-3     2022-04-03 [2] CRAN (R 4.2.1)
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
# Rsamtools            * 2.14.0    2022-11-01 [2] Bioconductor
# RSQLite                2.3.1     2023-04-03 [2] CRAN (R 4.2.3)
# rstudioapi             0.14      2022-08-22 [2] CRAN (R 4.2.1)
# rtracklayer            1.58.0    2022-11-01 [2] Bioconductor
# S4Vectors            * 0.36.2    2023-02-26 [2] Bioconductor
# scales                 1.2.1     2022-08-20 [2] CRAN (R 4.2.1)
# segmented              1.6-4     2023-04-13 [1] CRAN (R 4.2.3)
# sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.2.1)
# stringi                1.7.12    2023-01-11 [2] CRAN (R 4.2.2)
# stringr                1.5.0     2022-12-02 [2] CRAN (R 4.2.2)
# SummarizedExperiment * 1.28.0    2022-11-01 [2] Bioconductor
# survival               3.5-3     2023-02-12 [3] CRAN (R 4.2.3)
# tibble                 3.2.1     2023-03-20 [2] CRAN (R 4.2.3)
# tidyr                  1.3.0     2023-01-24 [2] CRAN (R 4.2.2)
# tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.2.1)
# tzdb                   0.3.0     2022-03-28 [2] CRAN (R 4.2.1)
# utf8                   1.2.3     2023-01-31 [2] CRAN (R 4.2.2)
# VariantAnnotation    * 1.44.1    2023-02-15 [2] Bioconductor
# vctrs                  0.6.2     2023-04-19 [1] CRAN (R 4.2.3)
# viridis              * 0.6.2     2021-10-13 [2] CRAN (R 4.2.1)
# viridisLite          * 0.4.2     2023-05-02 [1] CRAN (R 4.2.3)
# wesanderson          * 0.3.6     2018-04-20 [1] CRAN (R 4.2.2)
# WGCNA                * 1.72-1    2023-01-18 [1] CRAN (R 4.2.2)
# withr                  2.5.0     2022-03-03 [2] CRAN (R 4.2.1)
# xfun                   0.39      2023-04-20 [1] CRAN (R 4.2.3)
# XML                    3.99-0.14 2023-03-19 [2] CRAN (R 4.2.3)
# xml2                   1.3.3     2021-11-30 [2] CRAN (R 4.2.1)
# XVector              * 0.38.0    2022-11-01 [2] Bioconductor
# yaml                   2.3.7     2023-01-23 [2] CRAN (R 4.2.2)
# zlibbioc               1.44.0    2022-11-01 [2] Bioconductor
# 
# [1] /users/bsimbiat/R/4.2.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/library
# 
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────



