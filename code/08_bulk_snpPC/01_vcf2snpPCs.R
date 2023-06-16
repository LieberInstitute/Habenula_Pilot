#!/bin/env Rscript
## script to get MDS and snpPCs for a set of genotypes in a VCF.gz input data
## module plink must be loaded (or plink available in PATH)!
library(here)

vcf <- here("processed-data/08_bulk_snpPC/habenula_genotypes_n69.vcf.gz")

bedout <- sub(".vcf.gz", "", vcf, fixed = T)

## basic QC filter:
qcflt <- "--geno 0.1 --maf 0.05 --hwe 0.000001 --out"
## should not be applied for an already extracted VCF with samples of interest

## convert to plink BED:
cmd <- paste(
    "plink --make-bed --output-chr chrM --keep-allele-order --vcf",
    vcf, "--out", bedout
)
## or we can directly apply the QC filter during conversion
# cmd=paste("plink --make-bed --output-chr chrM --keep-allele-order --vcf",
#         vcf,  ,   bedout)
## for a file with 68 genotypes, this eliminates over half of the variants!
system(cmd)

indfile <- paste0(bedout, "_indep")

### independent and cluster
## --indep <window size>['kb'] <step size (variant ct)> <VIF threshold>
## produce a pruned subset of markers that are in approximate linkage equilibrium with each other
## --indep requires three parameters:
##          a window size in variant count or kilobase (if the 'kb' modifier is present) units,
##          a variant count to shift the window at the end of each step,
##          a variance inflation factor (VIF) threshold
cmd <- paste("plink --bfile ", bedout, "--indep 100 10 1.25 --out", bedout)
system(cmd)

## MDS components
# outmds=paste0(bedout, '_clmds')
cmd <- paste0(
    "plink --bfile ", bedout,
    " --cluster --mds-plot 10 --extract ", bedout, ".prune.in --out ", bedout
)
system(cmd)

# ## A transpose
cmd <- paste("plink --bfile", bedout, "--recode A-transpose --out", bedout)
system(cmd)

## read in genotypes
# genotypes  = read_delim(paste0(newbfile, ".traw"), delim="\t")
# snp = as.data.frame(genotypes[,-(1:6)])
# colnames(snp) = ifelse(grepl("^Br", ss(colnames(snp), "_")),
#                       ss(colnames(snp), "_"), ss(colnames(snp), "_",2))
# snp = as.matrix(snp[,unique(BrNums)])

#### read in MDS
mds <- read.table(paste0(bedout, ".mds"), header = TRUE, as.is = TRUE)

# rownames(mds) = ifelse(grepl("^Br", mds$FID),  mds$FID, mds$IID)
rmds <- mds[, -(1:3)] # remove FID, IID and SOL columns 1-3
colnames(rmds) <- paste0("snpPC", 1:ncol(rmds))
rmds$BrNum <- mds[,1]
data.table::setcolorder(rmds, "BrNum")

## write snpPCs file:
data.table::fwrite(rmds, file = paste0(bedout, ".snpPCs.tab"), sep = "\t")
