# November 4, 2022
# qc_plotter_Bukola.R - Using brainswapped rse objects for Habenula data
# to create QC plots for possible variable splitting.
# qrsh -l mem_free=50G,h_vmem=50G

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



## Loading data (brain swapped objects) ########################################
load(here("preprocessed_data", "count_data_bukola", 
          "rse_gene_Roche_Habenula_qcAndAnnotated_n69.Rdata")) # gene info
rse = rse_gene

# Checks:
table(rse$Flowcell, rse$PrimaryDx) 
table(rse$Sex, rse$PrimaryDx) 

### 1. Plotting "qc_plots_by.." ################################################
# For each QC metric plotted against diagnosis/flowcell and color coded by 
# diagnosis/flowcell.

### Base function: 
  create_boxplots <- function(objInt, cov_var, samp_cond, colorby){
    objInt = as.data.frame(colData(objInt))
    
    # creating df of possible titles
    orig_var_name <- c("ERCCsumLogErr", "numReads", "numMapped", "overallMapRate", 
                       "concordMapRate", "mitoRate", "totalAssignedGene", "rRNA_rate")
    var_plot_title <- c("ERCC RSS", "Num of Reads (log 10)", "Num Mapped (log 10)",
                        "Overall Map Rate", "Concordant Map Rate", "chrM Map Rate",
                        "Gene Assignment Rate", "Gene rRNA Rate")
    posib_title <- data.frame(orig_var_name, var_plot_title)
    
      if(cov_var %in% posib_title$orig_var_name) {
        newTitle <- posib_title[posib_title$orig_var_name == cov_var, 2]
      } else{
        newTitle <- cov_var
      }
      
    # In order to make sure geom_jitter and geom_text_repel use the same coordinates 
    # for points (to prevent mislabeling).
    pos <- position_jitter(seed = 2)
    
      plot = ggplot(objInt, aes_(x = objInt[,cov_var], y = as.factor(objInt[,samp_cond]))) +
        geom_boxplot(outlier.shape = NA) + 
        geom_jitter(aes_(color = as.factor(objInt[,colorby])), position = pos) +
        geom_text_repel(aes(label = objInt[,"BrNum"], color = 
                              as.factor(objInt[,colorby])), position = pos) +
                              theme_bw(base_size = 10) + 
        theme(legend.position= "top", plot.margin=unit (c (1.5,2,1,2), 'cm'), 
              axis.text.x = element_text(vjust = 0.7), text = element_text(size=15),
              axis.title = element_text(size=15)) +
        labs(x = newTitle, y = samp_cond) +
        guides(color = guide_legend(title = colorby))
  }

### Creating list of covariate names of interest:
covVarInt <-  c("ERCCsumLogErr", "numReads", "numMapped", "overallMapRate", 
                "concordMapRate", "mitoRate", "totalAssignedGene", "rRNA_rate")

## Creating plots for QC metrics against Flowcell.
  for(i in covVarInt){
    # function(objInt, cov_var, samp_cond, colorby)
    namer <- paste("plotflow", i, sep = "_")
    assign(namer, create_boxplots(rse, i, "Flowcell", "PrimaryDx"))   
  }
  # Print and save
    pdf("preprocessed_data/qc_qlots_bukola/qc_plots_byFlowCell.pdf", height = 7, width = 11)
    mget(ls(patt = "plotflow_"))
    dev.off()
    
## Creating plots for QC metrics against PrimaryDx.
  for(i in covVarInt){
    # function(objInt, cov_var, samp_cond, colorby)
    namer <- paste("plotdx", i, sep = "_")
    assign(namer, create_boxplots(rse, i, "PrimaryDx", "Flowcell"))   
  }
  # Print and save
    pdf("preprocessed_data/qc_qlots_bukola/qc_plots_byPrimaryDx.pdf", height = 7, width = 11)
    mget(ls(patt = "plotdx_"))
    dev.off()

#################################
# Based on smokingMouse pipeline:
#################################   
## Stable Variables ############################################################
# For ease, all downstream functions will utilize these variables.
  
  # pd & pd_dropped ####
pd = colData(rse)
pd = as.data.frame(pd)
  # dropped Race and Sex because all male and all Cauc.
drop = c("Brain.Region", "FQCbasicStats", "perBaseQual", "perTileQual",
           "GCcontent", "Ncontent", "SeqLengthDist", "SeqDuplication",
           "OverrepSeqs", "AdapterContent", "KmerContent", "SeqLength_R1",
           "perSeqQual", "perBaseContent", names(pd[,grepl("phred", names(pd))]),
           names(pd[,grepl("Adapter", names(pd))]), "SeqLength_R2", "bamFile",
           "trimmed", names(pd[,grepl("gene_", names(pd))]), "hasGenotype",
           "Age", "Race")
pd = pd[,!(names(pd)) %in% drop]
pd_dropped = colData(rse)[,(names(colData(rse))) %in% drop]

  # Log10 values
pd$numReads <- log10(pd$numReads)
pd$numMapped <- log10(pd$numMapped)
pd$numUnmapped <- log10(pd$numUnmapped)

  # QCmetCols ####
QCmetCols = c("RIN", "percentGC_R1", "percentGC_R2", "ERCCsumLogErr",
              "numReads", "numMapped", "numUnmapped", "overallMapRate",
              "concordMapRate", "totalMapped", "mitoMapped", "mitoRate",
              "totalAssignedGene", "rRNA_rate")

  # phenoCols ####
  # Creating intervals for AgeDeath
pd$AgeInterval = NA
levels = quantile(pd$AgeDeath, probs = c(0, 0.25, 0.5, 0.75, 1))
  
  for (i in 1:length(pd$AgeDeath)){
    if(levels[1] <= pd$AgeDeath[i] && pd$AgeDeath[i] < levels[2]){
      pd[i, "AgeInterval"] <- 1
    } else if(levels[2] <= pd$AgeDeath[i] && pd$AgeDeath[i] < levels[3]){
      pd[i, "AgeInterval"] <- 2
    } else if(levels[3] <= pd$AgeDeath[i] && pd$AgeDeath[i] < levels[4]){
      pd[i, "AgeInterval"] <- 3
    } else if(levels[4] <= pd$AgeDeath[i] && pd$AgeDeath[i] <= levels[5]){
      pd[i, "AgeInterval"] <- 4
    }}

phenoCols = c("AgeInterval", "PrimaryDx", "Flowcell")


  # rename_vars ####
# Creating df for plot text to rename variables:
orig_var_name <- c("RNum", "RIN", "BrNum", "AgeDeath", "Sex", "PrimaryDx", 
                   "percentGC_R1", "percentGC_R2", "ERCCsumLogErr", 
                   "numReads", "numMapped", "numUnmapped", "overallMapRate", 
                   "concordMapRate", "totalMapped", "mitoMapped", "mitoRate", 
                   "totalAssignedGene", "rRNA_rate", "Flowcell", "AgeInterval" )

var_plot_title <- c("RNum", "RIN", "Brain Number", "Age oof Death", "Sex",
                    "Primary Dx", "Percent GC R1", "Percent GC R2", 
                    "ERCC RSS", "Num of Reads (log 10)", "Num Mapped (log 10)",
                    "Num Unmapped (log10)", "Overall Map Rate", 
                    "Concordant Map Rate", "Total Mapped", "chrM Mapped",
                    "chrM Map Rate", "Total Assigned Genes", "Gene rRNA Rate",
                    "Flowcell", "Age Intervals")

rename_vars <- data.frame(orig_var_name, var_plot_title)

### 2. Plotting "Mito_vs_Ribo_" ################################################
# For each QC metric plotted against diagnosis/flowcell and color coded by 
# diagnosis/flowcell.
for (i in 1:length(phenoCols)){
  print(i)
  pheno_var = phenoCols[i]
  print(phenoCols[i])
  namer = paste("mito_vs_ribo", pheno_var, sep = "_by")
  assign(namer, 
        ggplot(pd, aes(x = 100*(mitoRate), y = log10(rRNA_rate), 
        color = as.factor(pd[,pheno_var]))) + geom_point() +
        labs(x = "Ribosomal Counts", y = "Percentage of MT Counts") +
        guides(color = guide_legend(title = 
          rename_vars[rename_vars$orig_var_name == pheno_var,]$var_plot_title))
        )
  rm(pheno_var)
}

# Plot colors are stuck on FlowCell. ****
pdf("preprocessed_data/qc_qlots_bukola/Mito_vs_Ribo_byPhenotype.pdf")
  grid.arrange(mito_vs_ribo_byAgeInterval, mito_vs_ribo_byFlowcell,
               mito_vs_ribo_byPrimaryDx, ncol = 1, 
               top=textGrob("Mito Rate vs Ribo Rate by Phenotype"))
dev.off()



    
    
    
## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
    
    

