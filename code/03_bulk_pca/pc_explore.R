##### 
# pca on many regions 
#####
library(SummarizedExperiment)
library(recount)
library(jaffelab)
library(RColorBrewer)

dir.create("pca_plots")
set.seed(5235)

## load joint data that has habenula
load("/dcl01/ajaffe/data/lab/brain_swap/count_data/RNAseq_Collection_postQC_n5544_11dataset_2020-09-07_geneRSE.Rdata")

## filter to relevant studies and samples
rel_dataset = c("BrainSeq_Phase2_DLPFC",
	"BrainSeq_Phase2_HIPPO", "BrainSeq_Phase3_Caudate", "BrainSeq_Phase4and5",
	"Habenula", "Nicotine_NAc", "psychENCODE_BP", "psychENCODE_MDD",
	"VA_PTSD")
keepIndex = which(rse$Dataset %in% rel_dataset & 
	rse$Sex == "M" & rse$Age > 20 & rse$Age < 69)
rse_subset = rse[,keepIndex]

## subset to same brains 
brain_tab = table(rse_subset$BrNum)
rse_hab = rse_subset[,rse_subset$Dataset == "Habenula"]
table(brain_tab[rse_hab$BrNum])
rse_subset_same = rse_subset[,rse_subset$BrNum %in% rse_hab$BrNum]

## add more samples to get to equal numbers
tt = table(rse_subset_same$Dataset)
rse_subset_notsame = rse_subset[, ! rse_subset$BrNum %in% rse_hab$BrNum]

## sample
num_to_sample = ncol(rse_hab) - tt
more_sample_ind = unlist(lapply(seq(along=num_to_sample), function(i) {
	sample(which(rse_subset_notsame$Dataset == names(num_to_sample)[i]), 
		num_to_sample[i], replace=FALSE)
}))

## recombine
rse_subset_even = cbind(rse_subset_same, rse_subset_notsame[,more_sample_ind])

## subset for expression
geneRpkm = getRPKM(rse_subset_even, "Length")
geneExprs_filter = log2(geneRpkm[rowMeans(geneRpkm) > 0.1,] +1 )

pca=  prcomp(t(geneExprs_filter))
pcaVars = getPcaVars(pca)

## plots
pdf("pca_plots/pc1_vs_pc2.pdf")
palette("paired" 