###

library(jaffelab)
library(SummarizedExperiment)
library(VariantAnnotation)

## load pipeline output
load("preprocessed_data/rse_gene_Roche_Habenula_PairedEnd_n73.Rdata")
load("preprocessed_data/rse_exon_Roche_Habenula_PairedEnd_n73.Rdata")
load("preprocessed_data/rse_jx_Roche_Habenula_PairedEnd_n73.Rdata")
load("preprocessed_data/rse_tx_Roche_Habenula_PairedEnd_n73.Rdata")

## make folders
dir.create("qc_qlots")
dir.create("count_data")

## original flow cell
pd = colData(rse_gene)
man = read.delim("preprocessed_data/.samples_unmerged.manifest",as.is=TRUE,header=FALSE)
man$Flowcell = ss(ss(man$V1, "/",8), "_",2)
pd$Flowcell = man$Flowcell[match(pd$SAMPLE_ID, man$V5)]

## technical checks
pdf("qc_qlots/technical_covariates_by_flowcell.pdf")
par(mar=c(5,6,2,2), cex.axis=1.8,cex.lab=1.8)
boxplot(pd$ERCCsumLogErr ~ pd$Flowcell,ylab="ERCC RSS",xlab="Flowcell")
boxplot(log10(pd$numReads)~ pd$Flowcell,ylab="log10(#Reads)",xlab="Flowcell")
boxplot(log10(pd$numMapped)~ pd$Flowcell,ylab="log10(#Align",xlab="Flowcell")
boxplot(pd$overallMapRate ~ pd$Flowcell,ylab="Overall Map Rate",xlab="Flowcell")
boxplot(pd$concordMapRate ~ pd$Flowcell,ylab="Concordant Map Rate",xlab="Flowcell")
boxplot(pd$mitoRate ~ pd$Flowcell,ylab="chrM Map Rate",xlab="Flowcell")
boxplot(pd$totalAssignedGene ~ pd$Flowcell,ylab="Gene Assignment Rate",xlab="Flowcell")
boxplot(pd$rRNA_rate ~ pd$Flowcell,ylab="Gene rRNA Rate",xlab="Flowcell")
dev.off()

# pheno data
pheno = read.csv("habenula_pheno_data.csv", as.is=TRUE)
pd = cbind(pheno[match(pd$SAMPLE_ID, pheno$RNum),], pd[,-1])
pd$BrNum[pd$BrNum == "Br0983"] = "Br983"

## more checks
table(pd$Flowcell, pd$PrimaryDx) # balanced
table(pd$Sex, pd$PrimaryDx) # all Male

## technical checks by dx
pdf("qc_qlots/technical_covariates_by_dx.pdf")
par(mar=c(5,6,2,2), cex.axis=1.8,cex.lab=1.8)
boxplot(pd$ERCCsumLogErr ~ pd$PrimaryDx,ylab="ERCC RSS",xlab="")
boxplot(log10(pd$numReads)~ pd$PrimaryDx,ylab="log10(#Reads)",xlab="")
boxplot(log10(pd$numMapped)~ pd$PrimaryDx,ylab="log10(#Align)",xlab="")
boxplot(pd$overallMapRate ~ pd$PrimaryDx,ylab="Overall Map Rate",xlab="")
boxplot(pd$concordMapRate ~ pd$PrimaryDx,ylab="Concordant Map Rate",xlab="")
boxplot(pd$mitoRate ~ pd$PrimaryDx,ylab="chrM Map Rate",xlab="")
boxplot(pd$totalAssignedGene ~ pd$PrimaryDx,ylab="Gene Assignment Rate",xlab="")
boxplot(pd$rRNA_rate ~ pd$PrimaryDx,ylab="Gene rRNA Rate",xlab="")
boxplot(pd$RIN ~ pd$PrimaryDx,ylab="RIN",xlab="")
boxplot(pd$AgeDeath ~ pd$PrimaryDx,ylab="Age",xlab="")
dev.off()

####################
## Genotype data ###
####################

## read in VCF
vcf = readVcf("preprocessed_data/Genotypes/mergedVariants.vcf.gz", "hg38" )
colnames(vcf) = ss(ss(colnames(vcf), "/", 8), "_")
vcf = vcf[,pd$RNum] # match up to correct BAM file

# add rs num
snpMap = import("/dcl01/ajaffe/data/lab/brain_swap/common_missense_SNVs_hg38.bed")
oo = findOverlaps(vcf, snpMap, type="equal")
rowData(vcf)$snpRsNum = NA
rowData(vcf)$snpRsNum[queryHits(oo)] = snpMap$name[subjectHits(oo)]
rowData(vcf)$snpRsNum[is.na(rowData(vcf)$snpRsNum)] = rownames(vcf)[is.na(rowData(vcf)$snpRsNum)]
rownames(vcf) = rowData(vcf)$snpRsNum

## filter
vcf = vcf[info(vcf)$DP > 5*ncol(vcf) & info(vcf)$DP < 80 *ncol(vcf) &
          nchar(ref(vcf)) == 1 & elementNROWS(alt(vcf)) == 1 &
          info(vcf)$VDB >0.1,]
		  
## read in obs
genotyped = readVcf("/dcl01/ajaffe/data/lab/brain_swap/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_QuadsPlus_GenotypingBarcode.vcf.gz")
br = ss(colnames(genotyped),"_")
pd$hasGenotype = pd$BrNum %in% br
table(pd$hasGenotype)

pd[!pd$BrNum %in% br,c("BrNum", "RNum")]
## add genotype info

mm = match(rownames(vcf), rownames(genotyped))
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

# > lapply(pdList, function(x) x[order(x$totalAssignedGene,decreasing=TRUE)[1:5],c(1:8, 49:54)])
# $Control
# DataFrame with 5 rows and 14 columns
         # RNum       RIN Brain.Region       BrNum         AgeDeath         Sex
  # <character> <numeric>  <character> <character>        <numeric> <character>
# 1      R18408       7.7     Habenula      Br5558            24.21           M
# 2      R18415       7.8     Habenula      Br5871 60.0547570157426           M
# 3      R18394       7.6     Habenula      Br2015            45.09           M
# 4      R18418       7.5     Habenula      Br8050 43.7151494410221           M
# 5      R18405       7.5     Habenula      Br5385            66.03           M
         # Race   PrimaryDx overallMapRate concordMapRate totalMapped mitoMapped
  # <character> <character>      <numeric>      <numeric>   <numeric>  <numeric>
# 1        CAUC     Control         0.7958         0.7617    82867834    3861781
# 2        CAUC     Control         0.9126         0.8708   333552410   12475778
# 3        CAUC     Control         0.8853         0.8366   336819020   11307680
# 4        CAUC     Control         0.8997         0.8653   101713405    4167531
# 5        CAUC     Control         0.8883         0.8298    56925849    1407553
            # mitoRate totalAssignedGene
           # <numeric>         <numeric>
# 1 0.0445266706187961 0.610463430908278
# 2 0.0360542245766406 0.580814794729586
# 3 0.0324815074511665 0.572761759870304
# 4 0.0393605417315162 0.570697391205041
# 5 0.0241294515961884 0.566993720739398

# $Schizo
# DataFrame with 5 rows and 14 columns
         # RNum       RIN Brain.Region       BrNum  AgeDeath         Sex
  # <character> <numeric>  <character> <character> <numeric> <character>
# 1      R18383         8     Habenula      Br5446     23.14           M
# 2      R18347       7.5     Habenula      Br1016     20.22           M
# 3      R18357       7.6     Habenula      Br1383     58.96           M
# 4      R18376       7.6     Habenula      Br5488     33.36           M
# 5      R18358       7.4     Habenula      Br1427     66.81           M
         # Race   PrimaryDx overallMapRate concordMapRate totalMapped mitoMapped
  # <character> <character>      <numeric>      <numeric>   <numeric>  <numeric>
# 1        CAUC      Schizo         0.8961         0.8468   250611560    7433588
# 2        CAUC      Schizo         0.9036         0.8646    65614599    2788549
# 3        CAUC      Schizo         0.8996         0.8492   301320099    8858087
# 4        CAUC      Schizo         0.9011          0.842    52911619    1834699
# 5        CAUC      Schizo         0.8434         0.8033   188932004    7227570
            # mitoRate totalAssignedGene
           # <numeric>         <numeric>
# 1 0.0288073155322417 0.597954206098505
# 2 0.0407663840266533 0.575220231395459
# 3 0.0285580592053627 0.550160125169923
# 4 0.0335127377881376 0.543869134924215
# 5 0.0368453593807254 0.540923339751301

lapply(pdList, function(x) x$BrNum[order(x$totalAssignedGene,decreasing=TRUE)[1:5]])


##############################
## cut genotype data #########
##############################
