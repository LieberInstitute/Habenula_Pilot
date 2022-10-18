# October 5, 2022
# Using Erik's qc_data.R file and Daianna's 02_QC.R file, I am updating the 
# code for checking the techincal effects of our Habenula bulk RNAseq data. 
# qrsh -l mem_free=100G,h_vmem=100G

library(jaffelab)
library(SummarizedExperiment)
library(VariantAnnotation)
library(here)
library(ggplot2)
library(ggrepel)

# Loading data (pipeline output). 
load(here("preprocessed_data", "rse_gene_Roche_Habenula_PairedEnd_n73.Rdata")) #gene info
load(here("preprocessed_data", "rse_exon_Roche_Habenula_PairedEnd_n73.Rdata")) #exon data
load(here("preprocessed_data","rse_jx_Roche_Habenula_PairedEnd_n73.Rdata")) #junction data
load(here("preprocessed_data", "rse_tx_Roche_Habenula_PairedEnd_n73.Rdata")) #transcript data

# Makes folders.
dir.create(here("preprocessed_data","qc_qlots_bukola"))
dir.create(here("preprocessed_data","count_data_bukola"))

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

# More checks: All male and about the same amount of SCZ and control between flowcells.
table(pd$Flowcell, pd$PrimaryDx) 
table(pd$Sex, pd$PrimaryDx) 

# Fixing columns with log to better view trends
pd$numReads = log10(pd$numReads)
pd$numMapped = log10(pd$numMapped)

## Base function creation for creating boxplots by different variables.
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

  # In order to make sure geom_jitter and geom_text_repel use the same coordinates 
  # for points (to prevent mislabeling).
    pos <- position_jitter(seed = 2)
  
  plot = ggplot(objInt, aes_(x = objInt[,cov_var], y = as.factor(objInt[,samp_cond]))) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(aes_(color = as.factor(objInt[,colorby])), position = pos) +
    geom_text_repel(aes(label = objInt[,"BrNum"], color = as.factor(objInt[,colorby])), position = pos) +
    theme_bw(base_size = 10) + 
    theme(legend.position= "top", plot.margin=unit (c (1.5,2,1,2), 'cm'), 
          axis.text.x = element_text(vjust = 0.7), text = element_text(size=15),
          axis.title = element_text(size=15)) +
    labs(x = newTitle, y = samp_cond) +
    guides(color = guide_legend(title = colorby))
    
}

# creating list of covariate names of interet
covVarInt <-  c("ERCCsumLogErr", "numReads", "numMapped", "overallMapRate", 
                "concordMapRate", "mitoRate", "totalAssignedGene", "rRNA_rate")

# Covariates by Flowcell
for(i in covVarInt){
  # function(objInt, cov_var, samp_cond, colorby)
  namer <- paste("plotflow", i, sep = "_")
  assign(namer, create_boxplots(pd, i, "Flowcell", "PrimaryDx"))   
}

  # printing plots
pdf("preprocessed_data/qc_qlots_bukola/qc_qlots_byFlowCell.pdf", height = 7, width = 11)
  mget(ls(patt = "plotflow_"))
dev.off()


# Covariates by Primary Diagnosis
for(i in covVarInt){
  # function(objInt, cov_var, samp_cond, colorby)
  namer <- paste("plotdx", i, sep = "_")
  assign(namer, create_boxplots(pd, i, "PrimaryDx", "Flowcell"))   
}

# printing plots
pdf("preprocessed_data/qc_qlots_bukola/qc_qlots_byPrimaryDx.pdf", height = 7, width = 11)
  mget(ls(patt = "plotdx_"))
dev.off()

####################
## Genotype data ###
####################

# Read in VCF
vcf = readVcf("preprocessed_data/Genotypes/mergedVariants.vcf.gz", "hg38" )
colnames(vcf) = ss(ss(colnames(vcf), "/", 8), "_")
