#!/bin/env Rscript
## script to get MDS and snpPCs for a set of genotypes in a VCF.gz input data
## module plink must be loaded (or plink available in PATH)!
library("here")
library("sessioninfo")



################################ Set up paths #################################

vcf <- '/dcs05/lieber/liebercentral/libdGenotype_LIBD001/BrainGenotyping/subsets/Habenula_n69/Hb_gt_merged_R.9_MAF.05_ann.vcf.gz'
out_dir = here('processed-data', '08_bulk_snpPC')

bedout <- file.path(out_dir, sub(".vcf.gz", "", basename(vcf), fixed = T))
bedout_filt <- file.path(
    out_dir, sub(".vcf.gz", "_filt", basename(vcf), fixed = T)
)

###############################################################################



############################ Convert to plink BED #############################

cmd <- paste(
    "plink --make-bed --output-chr chrM --keep-allele-order --geno 0.05 --hwe 1e-6 --maf 0.05 --vcf",
    vcf, "--out", bedout
)
system(cmd)

###############################################################################



############################ Delete sample Br5572 #############################

cmd <- paste(
    "plink --bfile", bedout,
    "--remove samples_to_drop.txt --make-bed --out", bedout_filt
)
system(cmd)

###############################################################################


########################### Independent and cluster ###########################

## --indep <window size>['kb'] <step size (variant ct)> <VIF threshold>
## produce a pruned subset of markers that are in approximate linkage equilibrium with each other
## --indep requires three parameters:
##          a window size in variant count or kilobase (if the 'kb' modifier is present) units,
##          a variant count to shift the window at the end of each step,
##          a variance inflation factor (VIF) threshold

cmd <- paste("plink --bfile ", bedout_filt, "--indep 100 10 1.25 --out", bedout_filt)
system(cmd)

###############################################################################



################################ MDS components ###############################

# outmds=paste0(bedout, '_clmds')
cmd <- paste0(
    "plink --bfile ", bedout_filt,
    " --cluster --mds-plot 10 --extract ", bedout_filt, ".prune.in --out ", bedout_filt
)
system(cmd)

# ## A transpose
cmd <- paste("plink --bfile", bedout_filt, "--recode A-transpose --out", bedout_filt)
system(cmd)

###############################################################################



############################### Extract SNP PCs ###############################

# genotypes  = read_delim(paste0(newbfile, ".traw"), delim="\t")
# snp = as.data.frame(genotypes[,-(1:6)])
# colnames(snp) = ifelse(grepl("^Br", ss(colnames(snp), "_")),
#                       ss(colnames(snp), "_"), ss(colnames(snp), "_",2))
# snp = as.matrix(snp[,unique(BrNums)])

#### read in MDS
mds <- read.table(paste0(bedout_filt, ".mds"), header = TRUE, as.is = TRUE)

# rownames(mds) = ifelse(grepl("^Br", mds$FID),  mds$FID, mds$IID)
rmds <- mds[, -(1:3)] # remove FID, IID and SOL columns 1-3
colnames(rmds) <- paste0("snpPC", 1:ncol(rmds))
rmds$BrNum <- mds[, 1]
data.table::setcolorder(rmds, "BrNum")

## write snpPCs file:
data.table::fwrite(rmds, file = paste0(bedout_filt, ".snpPCs.tab"), sep = "\t")

###############################################################################



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
