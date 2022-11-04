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
library(cowplot)
library(ggrepel)
library(gridExtra)
library(sessioninfo)

# Loading data (brain swapped objects)
load(here("preprocessed_data", "count_data_bukola", 
          "rse_gene_Roche_Habenula_qcAndAnnotated_n69.Rdata")) # gene info
rse = rse_gene

# Checks:
table(rse$Flowcell, rse$PrimaryDx) 
table(rse$Sex, rse$PrimaryDx) 

## Base function creation for creating boxplots by different variables.
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
    geom_text_repel(aes(label = objInt[,"BrNum"], color = as.factor(objInt[,colorby])), position = pos) +
    theme_bw(base_size = 10) + 
    theme(legend.position= "top", plot.margin=unit (c (1.5,2,1,2), 'cm'), 
          axis.text.x = element_text(vjust = 0.7), text = element_text(size=15),
          axis.title = element_text(size=15)) +
    labs(x = newTitle, y = samp_cond) +
    guides(color = guide_legend(title = colorby))
}

# Creating list of Covariate names of interest.
covVarInt <-  c("ERCCsumLogErr", "numReads", "numMapped", "overallMapRate", 
                "concordMapRate", "mitoRate", "totalAssignedGene", "rRNA_rate")
##############################################

## Covariates by Flowcell
  for(i in covVarInt){
    # function(objInt, cov_var, samp_cond, colorby)
    namer <- paste("plotflow", i, sep = "_")
    assign(namer, create_boxplots(rse, i, "Flowcell", "PrimaryDx"))   
  }
    # printing plots
    pdf("preprocessed_data/qc_qlots_bukola/qc_plots_byFlowCell.pdf", height = 7, width = 11)
    mget(ls(patt = "plotflow_"))
    dev.off()
#########################
    
## Covariates by Primary Diagnosis
  for(i in covVarInt){
    # function(objInt, cov_var, samp_cond, colorby)
    namer <- paste("plotdx", i, sep = "_")
    assign(namer, create_boxplots(rse, i, "PrimaryDx", "Flowcell"))   
  }
    # printing plots
    pdf("preprocessed_data/qc_qlots_bukola/qc_plots_byPrimaryDx.pdf", height = 7, width = 11)
    mget(ls(patt = "plotdx_"))
    dev.off()

#########################
# based on smokingMouse pipeline:
    
# stable variables
pd = colData(rse)
drop = c("Brain.Region", "FQCbasicStats", "perBaseQual", "perTileQual",
           "GCcontent", "Ncontent", "SeqLengthDist", "SeqDuplication",
           "OverrepSeqs", "AdapterContent", "KmerContent", "SeqLength_R1",
           "perSeqQual", "perBaseContent", names(pd[,grepl("phred", names(pd))]),
           names(pd[,grepl("Adapter", names(pd))]), "SeqLength_R2", "bamFile",
           "trimmed", names(pd[,grepl("gene_", names(pd))]), "hasGenotype")

pd = pd[,!(names(pd)) %in% drop]
pd_dropped = colData(rse)[,(names(colData(rse))) %in% drop]


# Mito Rate vs Ribo Rate ("mitoRate" (change to perc) vs "rRNA_rate")
mito_vs_ribo = pd[,c("mitoRate", "rRNA_rate", "Flowcell", "PrimaryDx", 
                     "Sex", "Race", "AgeDeath", 
                     )]
for(){
  
}    
    
mito_vs_ribo<- function (pheno_var, tissue, age, label){
  if (is.null(age)){
      RSE<-eval(parse_expr(paste("rse_gene_", tissue, sep="")))
    }
    else {
      RSE<-eval(parse_expr(paste("rse_gene", tissue, age, sep="_")))
     }
      if (label==""){
        plot=ggplot(data=as.data.frame(colData(RSE)), 
                    aes(x=subsets_Mito_percent, y=subsets_Ribo_percent, color=eval(parse_expr(pheno_var)), 
                        label=label))+ 
          geom_point() +
          theme(text = element_text(size = 10)) +
          theme(legend.position="right", plot.margin=unit (c (1.5,2,1,2), 'cm')) +
          labs(x="Percentage of mt counts", y="Percentage of ribosomal counts", color=pheno_var)
      }
      else {
        plot=ggplot(data=as.data.frame(colData(RSE)), 
                    aes(x=subsets_Mito_percent, y=subsets_Ribo_percent, color=eval(parse_expr(pheno_var)), 
                        label=eval(parse_expr(label))))+ 
          geom_point() +
          ## Add samples' labels
          geom_text_repel(color="black", size=2, max.overlaps = 100) +
          theme(text = element_text(size = 10)) +
          theme(legend.position="right", plot.margin=unit (c (1.5,2,1,2), 'cm')) +
          labs(x="Percentage of mt counts", y="Percentage of ribosomal counts", color=pheno_var)
      }
      return(plot)
    }
    


    
    
    
## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
    
    

