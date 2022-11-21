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
         "Age", "Race", "subsets_Mito_sum.1", "subsets_Mito_detected.1",
         "subsets_Mito_percent.1", "subsets_Ribo_sum.1", "subsets_Ribo_detected.1",
         "subsets_Ribo_percent.1", "sum.1", "detected.1")

pd = pd[,!(names(pd)) %in% drop]

# Changing classes for plotting
pd$percentGC_R1 <- as.numeric(as.character(pd$percentGC_R1))
pd$percentGC_R2 <- as.numeric(as.character(pd$percentGC_R2))

# Log10 values
pd$logNumReads <- log10(pd$numReads)
pd$logNumMapped <- log10(pd$numMapped)
pd$logNumUnmapped <- log10(pd$numUnmapped)

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
              "numReads", "numMapped", "numUnmapped", "overallMapRate",
              "concordMapRate", "totalMapped", "mitoMapped", "mitoRate",
              "totalAssignedGene", "rRNA_rate", "sum", "detected",
              "subsets_Mito_sum", "subsets_Mito_detected", "subsets_Mito_percent",
              "subsets_Ribo_sum", "subsets_Ribo_detected", "subsets_Ribo_percent",
              "logNumReads", "logNumMapped", "logNumUnmapped")

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
  
  # grabbing p val
  pval = pairwise.t.test(pd[, qc_metter], pd[, pheno], p.adjust.method = "bonferroni")
  pval = signif(pval)
  
  # Use pos to ensure jitter and text_repel share coordinates (prevents mislabeling).  
  pos <- position_jitter(seed = 2)
  
  plot = ggplot(pd, aes_(y = pd[,qc_metter], x = as.factor(pd[,pheno]))) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(aes_(color = as.factor(pd[,colorby])), position = pos) +
    geom_text_repel(aes(label = pd[,"BrNum"], color = as.factor(pd[,colorby])),
                    position = pos) +
    theme_bw(base_size = 10) + 
    theme(legend.position= "top", plot.margin=unit (c (1.5,2,1,2), 'cm'), 
          axis.text.x = element_text(vjust = 0.7), text = element_text(size=15),
          axis.title = element_text(size=15)) +
    labs(x = pheno, y = titler) +
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
# For each QC metric plotted against diagnosis/flowcell and color coded by 
# diagnosis/Flowcell.

# Base Function
mito_vs_ribo <- function(phenos){
  
  ggplot(pd, aes(x = 100*(mitoRate), y = log10(rRNA_rate),
                 color = as.factor(pd[,phenos]))) + geom_point() +
    labs(x = "Ribosomal Counts", y = "Percentage of MT Counts") +
    guides(color = guide_legend(title = 
                                  rename_vars[rename_vars$orig_var_name == 
                                                phenos,]$var_plot_title))
}

# Plot
plot2 = lapply(phenoCols, mito_vs_ribo)

# Save
pdf(file = here("preprocessed_data", "qc_qlots_bukola", "Mito_vs_Ribo_byPhenotype.pdf"))
  plot_grid(plot2[[1]], plot2[[2]], plot2[[3]], ncol = 1, labels =
              "Mito vs Ribo Rates by Phenotype", rel_heights = c(.65,.35)) 
dev.off()

### 3. Plotting "boxplot_[QCmetric]_vs_[pheno]" ########################################
applyQC = QCmetCols[QCmetCols != "RIN"]
phenoCols = unlist(phenoCols)

boxplot_qc_pheno <- function(QC_mets, phenos){
  plottingpd = pd[, c(QC_mets, phenos)]
  plottingpd[, phenos] = as.factor(plottingpd[, phenos])
  
  plot = ggplot(plottingpd, aes_(y = plottingpd[, QC_mets], 
                x = plottingpd[, phenos], fill = plottingpd[, phenos])) + 
    geom_boxplot(color="red", fill="orange") + 
    xlab(rename_vars[rename_vars$orig_var_name == phenos,]$var_plot_title) +
    ylab(rename_vars[rename_vars$orig_var_name == QC_mets,]$var_plot_title) 
  
  return(plot)
}

# Plotting
for (i in applyQC){
  
  fileName = here("preprocessed_data", "qc_qlots_bukola", "02_boxplot_qc_by_pheno", 
                  paste("Boxplot", i, "vs_phenos.pdf", sep = "_"))
  
  plotter = lapply(phenoCols, FUN = boxplot_qc_pheno, QC_mets = i)
  
  pdf(file = fileName, width = 5, height =  10)
    do.call(grid.arrange, list(grobs = plotter, ncol = 1, 
       top = paste(rename_vars[rename_vars$orig_var_name == i,]$var_plot_title,
       "by Phenotype")))
  dev.off()
}

### 4. Scatterplotting "QC_vs_QC" #####################################################
# QCmetCols
# [1] "RIN"               "percentGC_R1"      "percentGC_R2"     
# [4] "ERCCsumLogErr"     "numReads"          "numMapped"        
# [7] "numUnmapped"       "overallMapRate"    "concordMapRate"   
# [10] "totalMapped"       "mitoMapped"        "mitoRate"         
# [13] "totalAssignedGene" "rRNA_rate"        

for (i in QCmetCols){
  qc1 = i
  fileName = here("preprocessed_data", "qc_qlots_bukola", paste(qc1, 
            "vs_QC_mets.pdf", sep = "_"))
  
  for (j in QCmetCols){
    qc2 = j
    scatterer = ggplot(pd, aes(qc1, qc2)) + geom_point() 
    pdf(file = fileName)
    do.call(grid.arrange, list(grobs = scatterer, ncol = 1, top = 
            paste(rename_vars[rename_vars$orig_var_name == i,]$var_plot_title,
            "vs Other QC Metrics")))
    
    dev.off()
  }
}


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
