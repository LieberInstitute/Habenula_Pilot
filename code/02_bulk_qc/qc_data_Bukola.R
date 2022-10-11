# October 5, 2022
# Using Erik's qc_data.R file and Daianna's 02_QC.R file, I am updating the 
# code for checking the techincal effects of our Habenula bulk RNAseq data. 
# qrsh -l mem_free=100G,h_vmem=100G

library(jaffelab)
library(SummarizedExperiment)
library(VariantAnnotation)
library(here)
library(ggplot2)

# Loading data (pipeline output). 
load(here("preprocessed_data", "rse_gene_Roche_Habenula_PairedEnd_n73.Rdata")) #gene info
load(here("preprocessed_data", "rse_exon_Roche_Habenula_PairedEnd_n73.Rdata")) #exon data
load(here("preprocessed_data","rse_jx_Roche_Habenula_PairedEnd_n73.Rdata")) #junction data
load(here("preprocessed_data", "rse_tx_Roche_Habenula_PairedEnd_n73.Rdata")) #transcript data

# Makes folders.
dir.create("qc_qlots_bukola")
dir.create("count_data_bukola")

# original flow cell [ask Louise about]
pd = colData(rse_gene)
man = read.delim(here("preprocessed_data", ".samples_unmerged.manifest"), 
                 as.is=TRUE,header=FALSE)
man$Flowcell = ss(ss(man$V1, "/",8), "_",2)
pd$Flowcell = man$Flowcell[match(pd$SAMPLE_ID, man$V5)]

# Phenotype Information
pheno = read.csv(here("preprocessed_data", "habenula_pheno_data.csv"), as.is=TRUE)
pd = cbind(pheno[match(pd$SAMPLE_ID, pheno$RNum),], pd[,-1])
# pd$BrNum[startsWith(pd$BrNum, "Br0") == TRUE]
pd$BrNum[pd$BrNum == "Br0983"] = "Br983" #Why is this important

# Code for boxplot creation of covariates by flowcell
  # look into other variables rather than flowcell 
create_boxplots_flowcell <- function(objInt, cov_var, yaxTit) {
 # I don't have age information so I can't filter by it. ***
  objInt = as.data.frame(objInt)
   plot = ggplot(objInt, aes(x = cov_var, y = Flowcell)) +
      geom_boxplot(outlier.shape = NA) + 
      geom_jitter() +
      theme_classic(base_size = 10) +
      theme(legend.position="none", plot.margin=unit (c (1.5,2,1,2), 'cm'), 
            axis.text.x = element_text(vjust = 0.45) ) +
      labs(x=yaxTit, y= "Flowcell") 
   
  return(plot)
}

## Base function creation for creating boxplots
create_boxplots <- function(objInt, cov_var, samp_cond, colorby){
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
  
  coloring <- 
  
  plot = ggplot(objInt, aes_(x = objInt[,cov_var], y = as.factor(objInt[,samp_cond]))) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter() +
    theme_bw(base_size = 10) + 
    theme(legend.position= "top", plot.margin=unit (c (1.5,2,1,2), 'cm'), 
          axis.text.x = element_text(vjust = 0.45) ) +
    labs(x = newTitle, y = samp_cond, color = colorby) 
    
}

# printing for base function
pdf("qc_qlots_bukola/tester.pdf")
  plot
dev.off()


# Creating boxplots and printing them all onto one pdf.
pdf("qc_qlots_bukola/technical_covariates_by_flowcell.pdf")
# par(mfrow = , cex.axis=1.8,cex.lab=1.8)
 create_boxplots_flowcell(pd, pd$ERCCsumLogErr, "ERCC RSS")
 create_boxplots_flowcell(pd, log10(pd$numReads), "Num of Reads (log 10)")
 create_boxplots_flowcell(pd, log10(pd$numMapped), "Num Mapped (log 10)")
 create_boxplots_flowcell(pd, pd$overallMapRate, "Overall Map Rate")
 create_boxplots_flowcell(pd, pd$concordMapRate, "Concordant Map Rate")
 create_boxplots_flowcell(pd, pd$mitoRate, "chrM Map Rate")
 create_boxplots_flowcell(pd, pd$totalAssignedGene, "Gene Assignment Rate")
 create_boxplots_flowcell(pd, pd$rRNA_rate, "Gene rRNA Rate")
dev.off()



# More checks: All male and about the same amount of SCZ and control between flowcells.
table(pd$Flowcell, pd$PrimaryDx) 
table(pd$Sex, pd$PrimaryDx) 

# Code for boxplot creation of covariates by disease
create_boxplots_dx <- function(objInt, cov_var, yaxTit) {
  # I don't have age information so I can't filter by it. ***
  objInt = as.data.frame(objInt)
  plot = ggplot(objInt, aes(x = cov_var, y = PrimaryDx)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter() +
    theme_bw(base_size = 10) +
    theme(legend.position="none", plot.margin=unit (c (1.5,2,1,2), 'cm'), 
          axis.text.x = element_text(vjust = 0.45) ) + 
    labs(x=yaxTit, y= "Condition") 
  
  return(plot)
}

# QC by disease
pdf("qc_qlots_bukola/technical_covariates_by_dx.pdf")
# par(mar=c(5,6,2,2), cex.axis=1.8,cex.lab=1.8)
  create_boxplots_dx(pd, pd$ERCCsumLogErr, "ERCC RSS")
  create_boxplots_dx(pd, log10(pd$numReads), "Reads (log 10)")
  create_boxplots_dx(pd, pd$numMapped, "# Aligns(no log 10)")
  create_boxplots_dx(pd, pd$overallMapRate, "Overall Map Rate")
  create_boxplots_dx(pd, pd$concordMapRate, "Concordant Map Rate")
  create_boxplots_dx(pd, pd$mitoRate, "chrM Map Rate")
  create_boxplots_dx(pd, pd$totalAssignedGene, "Gene Assignment Rate")
  create_boxplots_dx(pd, pd$rRNA_rate, "Gene rRNA Rate")
  create_boxplots_dx(pd, pd$RIN, "RIN")
  create_boxplots_dx(pd, pd$AgeDeath, "Age")
dev.off()


# QC by other features?
# Plot all points: outliers = FALSE, and add geom_point or geom_scatter
# rework function to create one base 
# use geom_bw instead of geom_classic for grid
# load and us gg repel 
# make graphs more visually interesting
# share in habenula chat
# plotting flowcell and make samples different colors based on diagnosis 
# plotting diagnoses and make samples look different by flowcell
# start reading elsewhere



