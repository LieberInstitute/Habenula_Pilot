# November 4, 2022
# qc_plotter_Bukola.R - Using brainswapped rse objects for Habenula data
# to create QC plots for possible variable splitting.
# qrsh -l mem_free=20G,h_vmem=20G

library(SummarizedExperiment)
library(recount)
library(edgeR)
library(here)
library(ggplot2)
library(rlang)
library(scater)
library(jaffelab)
library(ggrepel)
library(gridExtra)
library(sessioninfo)
library(dplyr)
library(gridExtra)
library(grid)
library(cowplot)
library(grid)
library(ggplotify)
library(rstatix)
library(ggpubr)

## Loading data (brain swapped and filtered) ###################################
# gene
load(here("processed-data", "02_bulk_qc", "count_data_bukola", 
          "rse_gene_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata"))


# # exon
# load(here("processed-data", "02_bulk_qc", "count_data_bukola", 
#           "rse_exon_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata")) 
# rse_exon = rse_exon_filt
# 
# # jx
# load(here("processed-data", "02_bulk_qc", "count_data_bukola", 
#           "rse_jx_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata")) 
# rse_jx = rse_jx_filt
# 
# # tx
# load(here("processed-data", "02_bulk_qc", "count_data_bukola", 
#           "rse_tx_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata"))
# rse_tx = rse_tx_filt


## Focusing on rse_gene ########################################################
# Checks:
table(rse_gene$Flowcell, rse_gene$PrimaryDx) 
table(rse_gene$Sex, rse_gene$PrimaryDx) 

## Creating Stable Variables ###################################################
# For ease, all downstream functions will utilize these variables.

# pd & pd_dropped ####
pd = as.data.frame(colData(rse_gene))

# Changing schizo to SCZD
pd$PrimaryDx <- recode_factor(pd$PrimaryDx, Schizo = "SCZD")

# dropped Race and Sex because all male and all Cauc.
drop = c("Brain.Region", "FQCbasicStats", "perBaseQual", "perTileQual",
         "GCcontent", "Ncontent", "SeqLengthDist", "SeqDuplication",
         "OverrepSeqs", "AdapterContent", "KmerContent", "SeqLength_R1",
         "perSeqQual", "perBaseContent", names(pd[,grepl("phred", names(pd))]),
         names(pd[,grepl("Adapter", names(pd))]), "SeqLength_R2", "bamFile",
         "trimmed", names(pd[,grepl("gene_", names(pd))]), "hasGenotype",
         "Age", "Race", names(pd[,grepl("subsets_", names(pd))]), "Sex", "sum")

pd = pd[,!(names(pd)) %in% drop]

# Changing classes for plotting
pd$percentGC_R1 <- as.numeric(as.character(pd$percentGC_R1))
pd$percentGC_R2 <- as.numeric(as.character(pd$percentGC_R2))

# Log10 values
pd$logNumReads <- log10(pd$numReads)
pd$logNumMapped <- log10(pd$numMapped)
pd$logNumUnmapped <- log10(pd$numUnmapped)

drop2 = c("numReads", "numMapped", "numUnmapped")
pd = pd[,!(names(pd)) %in% drop2]

# Creating intervals for AgeDeath
pd$AgeInterval = NA
levels = quantile(pd$AgeDeath, probs = c(0, 0.25, 0.5, 0.75, 1))

for (i in 1:length(pd$AgeDeath)){
  if(levels[1] <= pd$AgeDeath[i] && pd$AgeDeath[i] < levels[2]){
    pd[i, "AgeInterval"] <- "20 to 31"
  } else if(levels[2] <= pd$AgeDeath[i] && pd$AgeDeath[i] < levels[3]){
    pd[i, "AgeInterval"] <- "31.5 to 43.5"
  } else if(levels[3] <= pd$AgeDeath[i] && pd$AgeDeath[i] < levels[4]){
    pd[i, "AgeInterval"] <- "43.75 to 55.9"
  } else if(levels[4] <= pd$AgeDeath[i] && pd$AgeDeath[i] <= levels[5]){
    pd[i, "AgeInterval"] <- "60 to 68"
  }}

# Variable phenotypes in our data
phenoCols = as.vector(c("AgeInterval", "PrimaryDx", "Flowcell"))
## Making phenotypes of interest into factors
    for (i in phenoCols){
      pd[,i] <- as.factor(pd[,i])
    }

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

### 1. Plotting "qc_plots_by.." ################################################
# QC mets plots by diagnosis/flowcell and colored by vice versa.

# Base function: 
create_boxplots <- function(pd, qc_metter, pheno, colorby){
  
  titler = rename_vars[rename_vars$orig_var_name == qc_metter, "var_plot_title"]
  
  # grabbing p value
  if(length(levels(pd[, pheno])) > 2){
    
    pval = aov(pd[,qc_metter] ~ pd[,pheno], data = pd)
    pval = signif(unlist(summary(pval))["Pr(>F)1"])
      
  } else if(length(levels(pd[, pheno])) <= 2){
    
      pval = pairwise.t.test(pd[, qc_metter], pd[, pheno])
      pval = signif(pval$p.value)
  }
  
  # coloring p-values in red for significance 
  if (pval <= 0.05/length(QCmetCols)){
    sigColor = "red"
  } else {
    sigColor = "blue"
  }
  
  # Use pos to ensure jitter and text_repel share coordinates (prevents mislabeling).  
  pos <- position_jitter(seed = 2)
  
  plot = ggplot(pd, aes_(x = pd[,qc_metter], y = as.factor(pd[,pheno]))) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(aes_(color = as.factor(pd[,colorby])), position = pos) +
    geom_text_repel(aes(label = pd[,"BrNum"], color = as.factor(pd[,colorby])),
                    position = pos) +
    theme_bw(base_size = 10) + 
    labs(y = pheno, x = titler, caption = paste("p-value =", pval)) +
    theme(legend.position= "bottom", plot.margin=unit (c (1.5,2,1,2), 'cm'), 
          axis.text.x = element_text(vjust = 0.7), text = element_text(size=15),
          axis.title = element_text(size=15), 
          plot.caption = element_text(color = sigColor, face = "italic")) +
    guides(color = guide_legend(title = colorby)) 
  print(plot)
}

##################################


# Save for QC by PrimaryDx with Flow colors
pdf(here("plots","02_bulk_qc", "qc_plots_bukola", "boxplot_qc_by_pheno", "qc_by_dx_boxplot.pdf"), height = 7, width = 11)
  for (i in QCmetCols){
    create_boxplots(pd, i , "PrimaryDx", "Flowcell")
  }
dev.off()

# Save for QC by Flowcell with Dx colors
pdf(here("plots","02_bulk_qc", "qc_plots_bukola", "boxplot_qc_by_pheno", "qc_by_flow_boxplot.pdf"), height = 7, width = 11)
for (i in QCmetCols){
  create_boxplots(pd, i , "Flowcell", "PrimaryDx")
}
dev.off()

# Save for QC by AgeInt with Flow colors
pdf(here("plots","02_bulk_qc", "qc_plots_bukola", "boxplot_qc_by_pheno", "qc_by_age_flow_boxplot.pdf"), height = 7, width = 11)
for (i in QCmetCols){
  create_boxplots(pd, i , "AgeInterval", "Flowcell")
}
dev.off()

# Save for QC by AgeInt with Dx colors
pdf(here("plots","02_bulk_qc", "qc_plots_bukola", "boxplot_qc_by_pheno", "qc_by_age_dx_boxplot.pdf"), height = 7, width = 11)
for (i in QCmetCols){
  create_boxplots(pd, i , "AgeInterval", "PrimaryDx")
}
dev.off()

### 2. Plotting "Mito_vs_Ribo_byPhenotype" #####################################
# Using subsets data, plotting mito vs ribo rate by phenotype

# Base Function
mito_vs_ribo <- function(pheno){
  
  # Use pos to ensure jitter and text_repel share coordinates (prevents mislabeling).  
  pos <- position_jitter(seed = 4)
  
  ggplot(pd, aes(x = mitoRate, y = rRNA_rate)) + 
    geom_point() +
    geom_jitter(aes_(color = as.factor(pd[,pheno])), position = pos) +
    geom_text_repel(aes(label = pd[,"BrNum"], color = as.factor(pd[,pheno])),
                    position = pos) +
    labs(x = "Ribosomal Counts", y = "Percentage of MT Counts", title = 
           paste("Mito vs Ribo Rates by", rename_vars[rename_vars$orig_var_name == 
              pheno,]$var_plot_title, sep = " ")) +
    guides(color = guide_legend(title =  rename_vars[rename_vars$orig_var_name == 
                                pheno,]$var_plot_title))
}

# Plot
plot2 = lapply(phenoCols, FUN = mito_vs_ribo)

# Save
pdf(file = here("plots","02_bulk_qc", "qc_plots_bukola", "Mito_vs_Ribo_byPhenotype.pdf"), width = 10, height = 15)
  plot_grid(plot2[[1]], plot2[[2]], plot2[[3]], ncol = 1) 
dev.off()

### 3. Scatterplotting "QC_vs_QC" #####################################################
boxplot_qc_qc <- function(QC1, QC2, pheno){
  
  pos <- position_jitter(seed = 7)
  
  plot = ggplot(pd, aes_(y = pd[, QC1], x = pd[, QC2])) + 
    geom_point() + 
    geom_jitter(aes_(color = as.factor(pd[,pheno])), position = pos) +
    geom_text_repel(aes(label = pd[,"BrNum"], color = as.factor(pd[,pheno])),
                    position = pos) +
    labs(x = rename_vars[rename_vars$orig_var_name == QC2,]$var_plot_title, 
    y = rename_vars[rename_vars$orig_var_name == QC1,]$var_plot_title) +
    guides(color = guide_legend(title =  rename_vars[rename_vars$orig_var_name == 
          pheno,]$var_plot_title))
  
  return(plot)
}

# Plotting Overall Map Rate vs Total Assigned Gene
pdf(file = here("plots", "02_bulk_qc", "qc_plots_bukola", "MapRate_vs_GeneAssigned.pdf"), width = 10, height = 15)
  plotqcs1 = lapply(phenoCols, FUN = boxplot_qc_qc, QC1 = "overallMapRate", QC2 = "totalAssignedGene")
  plot_grid(plotqcs1[[1]], plotqcs1[[2]], plotqcs1[[3]], ncol = 1)
dev.off()

# Plotting Detected vs rRNA Rate
pdf(file = here("plots", "02_bulk_qc", "qc_plots_bukola", "Detected_vs_rRNArate.pdf"), width = 10, height = 15)
  plotqcs2 = lapply(phenoCols, FUN = boxplot_qc_qc, QC1 = "detected", QC2 = "rRNA_rate")
  plot_grid(plotqcs2[[1]], plotqcs2[[2]], plotqcs2[[3]], ncol = 1)
dev.off()

# Plotting log10s of Num Mapped and Num Reads
pdf(file = here("plots", "02_bulk_qc", "qc_plots_bukola", "NumMapped_vs_NumReads.pdf"), width = 10, height = 15)
  plotqcs3 = lapply(phenoCols, FUN = boxplot_qc_qc, QC1 = "logNumReads", QC2 = "logNumMapped")
  plot_grid(plotqcs3[[1]], plotqcs3[[2]], plotqcs3[[3]], ncol = 1)
dev.off()

# Plotting Total Mapped vs Mito Mapped
pdf(file = here("plots", "02_bulk_qc", "qc_plots_bukola", "TotalMapped_vs_MitoMapped.pdf"), width = 10, height = 15)
  plotqcs4 = lapply(phenoCols, FUN = boxplot_qc_qc, QC2 = "mitoMapped", QC1 = "totalMapped")
  plot_grid(plotqcs4[[1]], plotqcs4[[2]], plotqcs4[[3]], ncol = 1)
dev.off()

# Plotting Overall Map Rate vs Concordant Map Rate
pdf(file = here("plots", "02_bulk_qc", "qc_plots_bukola", "Overall_vs_ConcordMapRate.pdf"), width = 10, height = 15)
  plotqcs5 = lapply(phenoCols, FUN = boxplot_qc_qc, QC2 = "overallMapRate", QC1 = "concordMapRate")
  plot_grid(plotqcs5[[1]], plotqcs5[[2]], plotqcs5[[3]], ncol = 1)
dev.off()

# Plotting Overall Map Rate vs Concordant Map Rate
pdf(file = here("plots", "02_bulk_qc", "qc_plots_bukola", "Overall_vs_ConcordMapRate.pdf"), width = 10, height = 15)
  plotqcs5 = lapply(phenoCols, FUN = boxplot_qc_qc, QC2 = "overallMapRate", QC1 = "concordMapRate")
  plot_grid(plotqcs5[[1]], plotqcs5[[2]], plotqcs5[[3]], ncol = 1)
dev.off()

# Plotting ERCC vs Mito Mapped Rate
pdf(file = here("plots", "02_bulk_qc", "qc_plots_bukola", "ERCC_vs_MitoMapRate.pdf"), width = 10, height = 15)
  plotqcs6 = lapply(phenoCols, FUN = boxplot_qc_qc, QC2 = "ERCCsumLogErr", QC1 = "mitoRate")
  plot_grid(plotqcs6[[1]], plotqcs6[[2]], plotqcs6[[3]], ncol = 1)
dev.off()

# Plotting Total Mapped vs Mito Mapped
pdf(file = here("plots", "02_bulk_qc", "qc_plots_bukola", "Total_vs_MitoMapRate.pdf"), width = 10, height = 15)
  plotqcs7 = lapply(phenoCols, FUN = boxplot_qc_qc, QC2 = "totalMapped", QC1 = "mitoRate")
  plot_grid(plotqcs7[[1]], plotqcs7[[2]], plotqcs7[[3]], ncol = 1)
dev.off()

## Reproducibility information
print('Reproducibility information:')
options(width = 120)
session_info( )
# 
# [1] "Reproducibility information:"
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
# abind                  1.4-5     2016-07-21 [2] CRAN (R 4.2.1)
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
# Biostrings             2.66.0    2022-11-01 [2] Bioconductor
# bit                    4.0.5     2022-11-15 [2] CRAN (R 4.2.2)
# bit64                  4.0.5     2020-08-30 [2] CRAN (R 4.2.1)
# bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.2.1)
# blob                   1.2.4     2023-03-17 [2] CRAN (R 4.2.3)
# broom                  1.0.4     2023-03-11 [2] CRAN (R 4.2.3)
# BSgenome               1.66.3    2023-02-16 [2] Bioconductor
# bumphunter             1.40.0    2022-11-01 [2] Bioconductor
# cachem                 1.0.8     2023-05-01 [1] CRAN (R 4.2.3)
# car                    3.1-2     2023-03-30 [2] CRAN (R 4.2.3)
# carData                3.0-5     2022-01-06 [2] CRAN (R 4.2.1)
# checkmate              2.1.0     2022-04-21 [2] CRAN (R 4.2.1)
# cli                    3.6.1     2023-03-23 [2] CRAN (R 4.2.3)
# cluster                2.1.4     2022-08-22 [3] CRAN (R 4.2.3)
# codetools              0.2-19    2023-02-01 [3] CRAN (R 4.2.3)
# colorout             * 1.2-2     2023-02-09 [1] Github (jalvesaq/colorout@79931fd)
# colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.2.2)
# cowplot              * 1.1.1     2020-12-30 [2] CRAN (R 4.2.1)
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
# doRNG                  1.8.6     2023-01-16 [2] CRAN (R 4.2.2)
# downloader             0.4       2015-07-09 [2] CRAN (R 4.2.1)
# dplyr                * 1.1.1     2023-03-22 [2] CRAN (R 4.2.3)
# edgeR                * 3.40.2    2023-01-19 [2] Bioconductor
# evaluate               0.21      2023-05-05 [1] CRAN (R 4.2.3)
# fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.2.2)
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
# ggplotify            * 0.1.0     2021-09-02 [1] CRAN (R 4.2.2)
# ggpubr               * 0.6.0     2023-02-10 [2] CRAN (R 4.2.2)
# ggrepel              * 0.9.3     2023-02-03 [2] CRAN (R 4.2.2)
# ggsignif               0.6.4     2022-10-13 [2] CRAN (R 4.2.1)
# glue                   1.6.2     2022-02-24 [2] CRAN (R 4.2.1)
# googledrive            2.1.0     2023-03-22 [2] CRAN (R 4.2.3)
# gridExtra            * 2.3       2017-09-09 [2] CRAN (R 4.2.1)
# gridGraphics           0.5-1     2020-12-13 [1] CRAN (R 4.2.2)
# gtable                 0.3.3     2023-03-21 [2] CRAN (R 4.2.3)
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.2.1)
# Hmisc                  5.0-1     2023-03-08 [2] CRAN (R 4.2.3)
# hms                    1.1.3     2023-03-21 [2] CRAN (R 4.2.3)
# htmlTable              2.4.1     2022-07-07 [2] CRAN (R 4.2.1)
# htmltools              0.5.5     2023-03-23 [2] CRAN (R 4.2.3)
# htmlwidgets            1.6.2     2023-03-17 [2] CRAN (R 4.2.3)
# httr                   1.4.5     2023-02-24 [2] CRAN (R 4.2.2)
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
# rlang                * 1.1.1     2023-04-28 [1] CRAN (R 4.2.3)
# rmarkdown              2.21      2023-03-26 [2] CRAN (R 4.2.3)
# rngtools               1.5.2     2021-09-20 [2] CRAN (R 4.2.1)
# rpart                  4.1.19    2022-10-21 [3] CRAN (R 4.2.3)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.2.1)
# Rsamtools              2.14.0    2022-11-01 [2] Bioconductor
# RSQLite                2.3.1     2023-04-03 [2] CRAN (R 4.2.3)
# rstatix              * 0.7.2     2023-02-01 [2] CRAN (R 4.2.2)
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
# withr                  2.5.0     2022-03-03 [2] CRAN (R 4.2.1)
# xfun                   0.39      2023-04-20 [1] CRAN (R 4.2.3)
# XML                    3.99-0.14 2023-03-19 [2] CRAN (R 4.2.3)
# xml2                   1.3.3     2021-11-30 [2] CRAN (R 4.2.1)
# XVector                0.38.0    2022-11-01 [2] Bioconductor
# yaml                   2.3.7     2023-01-23 [2] CRAN (R 4.2.2)
# yulab.utils            0.0.6     2022-12-20 [1] CRAN (R 4.2.2)
# zlibbioc               1.44.0    2022-11-01 [2] Bioconductor
# 
# [1] /users/bsimbiat/R/4.2.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/library
# 
