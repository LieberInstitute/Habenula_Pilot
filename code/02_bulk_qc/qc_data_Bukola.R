# October 5, 2022
# Using Erik's qc_data.R file and Daianna's 02_QC.R file, I am updating the 
# code for checking the techincal effects of our Habenula bulk RNAseq data. 
# qrsh -l mem_free=100G,h_vmem=100GR

library(jaffelab)
library(SummarizedExperiment)
library(VariantAnnotation)
library(here)

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

# Code for boxplot creation
create_boxplots <- function(pheno_var, qc_var, tissue, age) {
  if (is.null(age)){
    ## Tissue data
    RSE<-eval(parse_expr(paste("rse_gene_", tissue, sep="")))
  }
  else {
    ## Tissue and Age data
    RSE<-eval(parse_expr(paste("rse_gene", tissue, age, sep="_")))
  }
  
  ## Quantitative QC values grouped by a qualitative phenotype variable
  plot=ggplot(as.data.frame(colData(RSE)), 
              aes(x=eval(parse_expr(pheno_var)), y=eval(parse_expr(qc_var)), 
                  fill=eval(parse_expr(pheno_var)))) +
    geom_boxplot() +
    theme_classic(base_size = 10) +
    theme(legend.position="none", plot.margin=unit (c (1.5,2,1,2), 'cm'), 
          axis.text.x = element_text(vjust = 0.45) ) +
    labs(x=pheno_var, y=qc_var) 
  return(plot)
}

# Creating boxplots and printing them all onto one pdf.
pdf("qc_qlots_bukola/technical_covariates_by_flowcell.pdf")
par(mar=c(5,6,2,2), cex.axis=1.8,cex.lab=1.8)
create_boxplots(ERCCsumLogErr, Flowcell)
  
)
dev.off()
