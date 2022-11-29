# November 29, 2022
# 01_pc_explore_Bukola.R - Creating  PCA plots for bulk RNA-seq data pre and post
# brain sample drops.
# qrsh -l mem_free=20G,h_vmem=20G

library(SummarizedExperiment)
library(recount)
library(jaffelab)
library(RColorBrewer)

## Adapted from Josh. Don't really like. 
rse = read.csv("/dcs04/lieber/lcolladotor/libdDataSwaps_LIBD001/brain_swap/RNAseq_Collection_postQC_n5544_11dataset_2020-09-07_phenotypes.csv")

## filter to relevant studies and samples
rel_dataset = c("BrainSeq_Phase2_DLPFC",
                "BrainSeq_Phase2_HIPPO", "BrainSeq_Phase3_Caudate", "BrainSeq_Phase4and5",
                "Habenula", "Nicotine_NAc", "psychENCODE_BP", "psychENCODE_MDD",
                "VA_PTSD")
keepIndex = which(rse$Dataset %in% rel_dataset & 
                    rse$Sex == "M" & rse$Age > 20 & rse$Age < 69)
rse_subset = rse[keepIndex,]

## subset to same brains 
brain_tab = table(rse_subset$BrNum)
rse_hab = rse_subset[rse_subset$Dataset == "Habenula",]
table(brain_tab[rse_hab$BrNum])
rse_subset_same = rse_subset[rse_subset$BrNum %in% rse_hab$BrNum,]

## Adapted from smokingMouse. Like much better:


## Testing pca generator:
pca<-prcomp(t(assays(rse_hab)$logcounts))

## % of the variance explained by each PC
pca_vars<- getPcaVars(pca)
pca_vars_labs<- paste0(
  "PC", seq(along = pca_vars), ": ",
  pca_vars, "% Var Expl")
