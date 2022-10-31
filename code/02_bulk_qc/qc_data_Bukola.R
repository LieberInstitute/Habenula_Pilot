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
# dir.create(here("preprocessed_data","qc_qlots_bukola"))
# dir.create(here("preprocessed_data","count_data_bukola"))

# original flow cell
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


### FROM GitHub: 

####################
## Genotype data ###
####################

# Read in VCF
vcf = readVcf("preprocessed_data/Genotypes/mergedVariants.vcf.gz", "hg38" )
colnames(vcf) = ss(basename(colnames(vcf)), "_")

# Match order of RNums 
all(pd$RNum %in% colnames(vcf))
# TRUE. Otherwise, find another method of subsetting
vcf = vcf[,pd$RNum] 

# add rs num
snpMap = import(here("/dcs04/lieber/lcolladotor/libdDataSwaps_LIBD001/brain_swap/common_missense_SNVs_hg38.bed"))
oo = findOverlaps(query = vcf, subject = snpMap, type="equal")
rowData(vcf)$snpRsNum = NA
rowData(vcf)$snpRsNum[queryHits(oo)] = snpMap$name[subjectHits(oo)]
any(is.na(rowData(vcf)$snpRsNum))
# FALSE. If not false, then we will need another method for rownames.
rownames(vcf) = rowData(vcf)$snpRsNum

## filter ****
vcf = vcf[info(vcf)$DP > 5*ncol(vcf) & info(vcf)$DP < 80 *ncol(vcf) &
            nchar(ref(vcf)) == 1 & elementNROWS(alt(vcf)) == 1 &
            info(vcf)$VDB >0.1,]

## read in obs
genotyped = readVcf("/dcs04/lieber/lcolladotor/libdDataSwaps_LIBD001/brain_swap/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_QuadsPlus_GenotypingBarcode.vcf.gz")
brainNumbers = ss(colnames(genotyped),"_")
pd$hasGenotype = pd$BrNum %in% brainNumbers
table(pd$hasGenotype)

# Finds brains without an ID number in the imputed data.
pd[!(pd$BrNum %in% brainNumbers), c("BrNum", "RNum")]


## add genotype info ****
matcher = match(rownames(vcf), rownames(genotyped))
vcf = vcf[!is.na(mm),]
genotyped = genotyped[mm[!is.na(mm)],]

## RNA-seq SNPs
snpsCalled = geno(vcf)$GT
snpsCalled[snpsCalled == "."] = 0
snpsCalled[snpsCalled == "0/1"] = 1
snpsCalled[snpsCalled == "1/1"] = 2
class(snpsCalled) = "numeric"		  

## DNA genotypes
snpsGeno = geno(genotyped)$GT
snpsGeno[snpsGeno == "."] = NA
snpsGeno[snpsGeno == "0/0"] = 0
snpsGeno[snpsGeno == "0/1"] = 1
snpsGeno[snpsGeno == "1/1"] = 2
class(snpsGeno) = "numeric"

## flip
table(rowRanges(vcf)$REF == rowRanges(genotyped)$REF  | 
        rowRanges(vcf)$REF == as.character(unlist(alt(genotyped))) )
toFlip = which(rowRanges(vcf)$REF != rowRanges(genotyped)$REF)
snpsGeno[toFlip,] = 2-snpsGeno[toFlip,] 

## correlate
colnames(snpsGeno)= ss(colnames(snpsGeno), "_")
snpCorObs = cor(snpsCalled, snpsGeno, use="pairwise.complete.obs")

## check best
bestCor = data.frame(maxIndex = apply(as.matrix(snpCorObs), 1, which.max),
                     maxCor = apply(snpCorObs, 1, max, na.rm=TRUE))
bestCor$RNA_BrNum = pd$BrNum[match(rownames(bestCor), pd$RNum)]
bestCor$DNA_BrNum = colnames(snpCorObs)[bestCor$maxIndex]
bestCor$maxIndex = NULL
bestCor$RNum = rownames(bestCor)
bestCor$RNA_Region = pd$Region[match(bestCor$RNum, pd$RNum)]

checkCor = bestCor[bestCor$RNA_BrNum != bestCor$DNA_BrNum & bestCor$maxCor > 0.6,]
checkCor

b = unique(c(checkCor$RNA_BrNum, checkCor$DNA_BrNum))
bestCor[bestCor$RNA_BrNum %in% b | bestCor$DNA_BrNum %in% b,]

bestCor[bestCor$DNA_BrNum == "Br1750",]

### drop some samples for identity/genotype
dropRNAs = c("R18355", "R18393", "R18421","R18364")
pd = pd[! pd$RNum %in% dropRNAs,]

## swap some samples
pd[pd$RNum == "R18422",4:8] = pheno[pheno$BrNum == "Br1350",4:8]
pd[pd$RNum == "R18423",4:8] = pheno[pheno$BrNum == "Br1225",4:8]

################################
############ save counts #######
################################

## save gene counts
rse_gene = rse_gene[,pd$RNum]
colData(rse_gene) = pd
save(rse_gene, file = paste0("count_data/rse_gene_Roche_Habenula_qcAndAnnotated_n",
                             ncol(rse_gene), ".Rdata"))

## save exons counts
rse_exon = rse_exon[,pd$RNum]
colData(rse_exon) = pd
save(rse_exon, file = paste0("count_data/rse_exon_Roche_Habenula_qcAndAnnotated_n",
                             ncol(rse_exon), ".Rdata"))

## save junction counts
rse_jx = rse_jx[,pd$RNum]
colData(rse_jx) = pd
## filter
jIndex = (rowSums(assays(rse_jx)$counts) > 9) & (rowData(rse_jx)$Class != "Novel")
rse_jx = rse_jx[jIndex,]
save(rse_jx, file = paste0("count_data/rse_jx_Roche_Habenula_qcAndAnnotated_n",
                           ncol(rse_jx), ".Rdata"))

## save tx counts
rse_tx = rse_tx[,pd$RNum]
colData(rse_tx) = pd
save(rse_tx, file = paste0("count_data/rse_tx_Roche_Habenula_qcAndAnnotated_n",
                           ncol(rse_tx), ".Rdata"))


#######################	
## rank on quality ####
pdList = split(pd, pd$PrimaryDx)
lapply(pdList, function(x) x[order(x$totalAssignedGene,decreasing=TRUE)[1:5],c(1:8, 49:54)])
