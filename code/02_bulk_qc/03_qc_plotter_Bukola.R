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
rse_gene = rse_gene_filt

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
Sys.time()
proc.time()
options(width = 120)
session_info()
