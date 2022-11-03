# October 5, 2022
# Using Erik's qc_data.R file and Daianna's 02_QC.R file, I am updating the 
# code for checking the techincal effects of our Habenula bulk RNAseq data. 
# qrsh -l mem_free=50G,h_vmem=50G

library(jaffelab)
library(SummarizedExperiment)
library(VariantAnnotation)
library(here)
library(ggplot2)
library(ggrepel)

# Loading data (pipeline output). 
load(here("preprocessed_data", "rse_gene_Roche_Habenula_PairedEnd_n73.Rdata")) # gene info
load(here("preprocessed_data", "rse_exon_Roche_Habenula_PairedEnd_n73.Rdata")) # exon data
load(here("preprocessed_data","rse_jx_Roche_Habenula_PairedEnd_n73.Rdata")) # junction data
load(here("preprocessed_data", "rse_tx_Roche_Habenula_PairedEnd_n73.Rdata")) # transcript data

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
pd$BrNum[startsWith(pd$BrNum, "Br0") == TRUE] 
pd$BrNum[pd$BrNum == "Br0983"] = "Br983" # Matching LIBD naming conventions.  

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

# Reads VCF. Contains info on SNPs imputed off of our bulk RNAseq data. 
vcf = readVcf("preprocessed_data/Genotypes/mergedVariants.vcf.gz", "hg38" )
colnames(vcf) = ss(basename(colnames(vcf)), "_") # Grabs name for each SNP from filename.

# Match order of RNums. Makes sure all samples in phenotyped data match the samples of our imputed data. 
all(pd$RNum %in% colnames(vcf))
# Reorders imputec samples to match the order of the phenotyped samples.
vcf = vcf[,pd$RNum] 

# Adding RsNum (SNP information). RsNum is a method of naming a specific SNP rather than calling on its chromosome and base location.
snpMap = import(here("/dcs04/lieber/lcolladotor/libdDataSwaps_LIBD001/brain_swap/common_missense_SNVs_hg38.bed"))
oo = findOverlaps(query = vcf, subject = snpMap, type="equal") # finds the matching RsNums between our hg38 snpMap and our imputed RNAseq data.
rowData(vcf)$snpRsNum = NA
rowData(vcf)$snpRsNum[queryHits(oo)] = snpMap$name[subjectHits(oo)] # Adds correct RsNum to each SNP for all brain samples.
any(is.na(rowData(vcf)$snpRsNum)) # checks to see if we've missed any SNP names.
# FALSE. If not false, then we will need another method for naming the vcf. 
rownames(vcf) = rowData(vcf)$snpRsNum

## filter
vcf = vcf[info(vcf)$DP > 5*ncol(vcf) & info(vcf)$DP < 80 *ncol(vcf) &
            nchar(ref(vcf)) == 1 & elementNROWS(alt(vcf)) == 1 &
            info(vcf)$VDB >0.1,]

## Load genotyped object that contains info on (brains that have already been genotyped across LIBD database).
genotyped = readVcf("/dcs04/lieber/lcolladotor/libdDataSwaps_LIBD001/brain_swap/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_QuadsPlus_GenotypingBarcode.vcf.gz")
brainNumbers = ss(colnames(genotyped),"_") # Grabs brain number for each sample from file name.
pd$hasGenotype = pd$BrNum %in% brainNumbers # Checks if our RNA seq samples can be found in the already genotyped database.
table(pd$hasGenotype) # 3 are not.

# Pulls the 3 and displays their brain number and sample name.
pd[pd$hasGenotype == FALSE, c("BrNum", "RNum")]

## Adding genotype information to  samples that lack it.
matcher = match(rownames(vcf), rownames(genotyped)) # Finds location of imputed SNP data among known genotyped object
vcf = vcf[!is.na(matcher),] # drops all brain samples that have SNP information that cannot be found in our genotyped pool
genotyped = genotyped[matcher[!is.na(matcher)],]  # subsets database of known genotypes for only samples that are present in our data 

## RNA-seq SNPs (bulk RNAseq data) imputed off of our RNA data
### Making a numerical system for each SNP for later correlative analyses.
snpsCalled = geno(vcf)$GT
snpsCalled[snpsCalled == "."] = 0 #homozygous reference
snpsCalled[snpsCalled == "0/1"] = 1 # heterozygous
snpsCalled[snpsCalled == "1/1"] = 2 # homozygous alternative
class(snpsCalled) = "numeric"		  

## DNA genotypes 
### Making a numerical system for each SNP for later correlative analyses.
### Idea: RNA SNPs should highly correlate to respective genotype which will confirm identity of unknown samples.
snpsGeno = geno(genotyped)$GT
snpsGeno[snpsGeno == "." | snpsGeno == "./."] = NA
snpsGeno[snpsGeno == "0/0"] = 0 # homozygous reference
snpsGeno[snpsGeno == "0/1"] = 1 # heterozygous
snpsGeno[snpsGeno == "1/1"] = 2 # homozygous alternative
class(snpsGeno) = "numeric"

## flip
## Checking to make sure that SNPs in VCF are either ref or alt bases. Any other value indicates major errors upstream.
table(rowRanges(vcf)$REF == rowRanges(genotyped)$REF  | 
        rowRanges(vcf)$REF == as.character(unlist(rowRanges(genotyped)$ALT))) # no FALSE is good.
toFlip = which(rowRanges(vcf)$REF == as.character(unlist(rowRanges(genotyped)$ALT)))
snpsGeno[toFlip,] = 2-snpsGeno[toFlip,] # Changing numerical values of imputed data where reference values match genotyped alt.

## Creating correlation matrix.
colnames(snpsGeno)= ss(colnames(snpsGeno), "_")
snpCorObs = cor(snpsCalled, snpsGeno, use="pairwise.complete.obs")

## check best
bestCor = data.frame(maxIndex = apply(as.matrix(snpCorObs), 1, which.max),
                     maxCor = apply(snpCorObs, 1, max, na.rm=TRUE)) # Grabs max correlation values from matrix as well as its index. 
bestCor$RNA_BrNum = pd$BrNum[match(rownames(bestCor), pd$RNum)] # Grabs brain number for each RNum found in imputed pheno data.  
bestCor$DNA_BrNum = colnames(snpCorObs)[bestCor$maxIndex] # Grabs brain number for each sample based on max correlation between sample and genotyped data.
bestCor$maxIndex = NULL
bestCor$RNum = rownames(bestCor)
# bestCor$RNA_Region = pd$Region[match(bestCor$RNum, pd$RNum)]  There is no pd$Region

# Subsets for area of possible brain swapping where we are confident of the correlation values.
checkCor = bestCor[bestCor$RNA_BrNum != bestCor$DNA_BrNum & bestCor$maxCor > 0.6,] 
checkCor
# Drop samples that are swapped but are not convincingly anything else.
dropRNAs = rownames(bestCor[bestCor$RNA_BrNum != bestCor$DNA_BrNum & bestCor$maxCor < 0.6,])
pd = pd[! pd$RNum %in% dropRNAs,]

# Checking the brain samples involved in mixup with high correlations to a known genotype.
b = unique(c(checkCor$RNA_BrNum, checkCor$DNA_BrNum)) 
bestCor[bestCor$RNA_BrNum %in% b | bestCor$DNA_BrNum %in% b,] 
#          maxCor RNA_BrNum DNA_BrNum   RNum
# R18420 0.8943435    Br5212    Br5212 R18420 # Makes sense. 
# R18421 0.9406773    Br1350    Br5212 R18421 # We already have this Br5212. No duplicates. [Delete R18421]
# R18422 0.9424066    Br1225    Br1350 R18422 
# R18423 0.9043557    Br1750    Br1225 R18423 

# Dropping sample R18421 
pd = pd[pd$RNum != "R18421",]

# Checking remaining samples
bestCor[bestCor$DNA_BrNum == "Br1750",] # This brain number doesn't exist in genotyped data.
bestCor[bestCor$DNA_BrNum == "Br1225",]
#           maxCor RNA_BrNum DNA_BrNum   RNum
# R18423 0.9043557    Br1750    Br1225 R18423 

## swap some samples
pd[pd$RNum == "R18422",4:8] = pheno[pheno$BrNum == "Br1350",4:8] # Br1225 <- Br1350
pd[pd$RNum == "R18423",4:8] = pheno[pheno$BrNum == "Br1225",4:8] # Br1750 <-Br1225

################################
############ save counts #######
################################

## save gene counts
rse_gene = rse_gene[,pd$RNum]
colData(rse_gene) = pd
save(rse_gene, file = here("preprocessed_data", paste0("count_data_bukola/rse_gene_Roche_Habenula_qcAndAnnotated_n",
                             ncol(rse_gene), ".Rdata")))

## save exons counts
rse_exon = rse_exon[,pd$RNum]
colData(rse_exon) = pd
save(rse_exon, file = here("preprocessed_data", paste0("count_data/rse_exon_Roche_Habenula_qcAndAnnotated_n",
                             ncol(rse_exon), ".Rdata")))

## save junction counts
rse_jx = rse_jx[,pd$RNum]
colData(rse_jx) = pd
## filter
jIndex = (rowSums(assays(rse_jx)$counts) > 9) & (rowData(rse_jx)$Class != "Novel")
rse_jx = rse_jx[jIndex,]
save(rse_jx, file = here("preprocessed_data", paste0("count_data/rse_jx_Roche_Habenula_qcAndAnnotated_n",
                           ncol(rse_jx), ".Rdata")))

## save tx counts
rse_tx = rse_tx[,pd$RNum]
colData(rse_tx) = pd
save(rse_tx, file = here("preprocessed_data", paste0("count_data/rse_tx_Roche_Habenula_qcAndAnnotated_n",
                           ncol(rse_tx), ".Rdata")))


#######################	
## rank on quality ####
pdList = split(pd, pd$PrimaryDx)
lapply(pdList, function(x) x[order(x$totalAssignedGene,decreasing=TRUE)[1:5],c(1:8, 49:54)])
