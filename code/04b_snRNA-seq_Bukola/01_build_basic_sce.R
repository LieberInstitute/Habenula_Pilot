# January 6, 2023 - Bukola Ajanaku
# Building basic sce objects. This follows the initial pca analysis for bulk 
# Habenula data. Only 7 samples, all control.
# Based on: 
# 1) https://github.com/LieberInstitute/Habenula_Bulk/blob/master/code/08_snRNA-seq_Erik/20210323_human_hb_neun.R
# 2) https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/code/03_build_sce/01_build_basic_sce.R
# 3) https://github.com/LieberInstitute/Habenula_Bulk/blob/master/code/09_snRNA-seq_re-processed/00_get_droplet_scores.R
# qrsh -l mem_free=50G,h_vmem=50G

# Loading relevant libraries
library("SingleCellExperiment")
library("SpatialExperiment")
library("DropletUtils")
library("here")
library("rtracklayer")
library("lobstr")
library("sessioninfo")
library("dplyr")
library("scater")

# Loading raw data:
load(here("processed-data","09_snRNA-seq_re-processed","20220601_human_hb_processing.rda"))

# Grabbing locations of raw data for reading in
addressRawData <- unique(sce.all.hb$Sample)

# Reading in data
sce_all = list()

for (i in 1:length(addressRawData)){
  sce_all[[i]] <- read10xCounts(addressRawData[i], col.names=TRUE) 
}

# Match rownames and appending by column
for(i in 2:length(sce_all)){
  sce_all[[i]] <- sce_all[[i]][match(rownames(sce_all[[1]]), rownames(sce_all[[i]])),]
}

sce_hb <- cbind(sce_all[[1]], sce_all[[2]], sce_all[[3]], sce_all[[4]], sce_all[[5]],
                sce_all[[6]], sce_all[[7]])
# 36601x7438664 dim. Much bigger than what Erik got (36601 x 3326836 dim)
rm(sce_all)

rownames(sce_hb) <- uniquifyFeatureNames(rowData(sce_hb)$ID, rowData(sce_hb)$Symbol)

## QUALITY CONTROL
# Thought process: (this is 10x Genomics data)
# 1) Based on OSCA book, can use mean absolute deviation (MADs) approach to 
# determine outliers but this is not as "straightforward" as it may through 
# away neurons.
# 2) There is a test using UMI/barcode rank (the knee plots), to determine 
# individual thresholds for what is too low of quality but this may drop cells 
# with organically low RNA content.
# 3) This is a single-nuc RNA project meaning that dropping high mito content
# is a useful QC metric because it shows samples where cytoplasm 
# was not fully or successfully stripped.
# Game plan:
# 1) Use subsets_Mito_percent to drop high mito content droplets.
# 2) Use DLPFC example of utilizing knee inflection plots to create further drop
# thresholds.


## 1) Discarding cells with high mito proportions:

# Saving unfiltered sce for future reference and plotting purposes.
unfiltered_sce <- sce_hb

# finding genes associated with mitochondrial genome
is.mito <- grep("^MT-", rowData(sce_hb)$Symbol)
stats <- perCellQCMetrics(sce_hb, subsets = list(Mito = is.mito))

# grabbing droplets with relatively high mito (based on MADs, no preset threshold)
high.mito <- isOutlier(stats$subsets_Mito_percent, type = "higher")

# dropping high mito content
discard <- high.mito
sce_hb <- sce_hb[,!discard]

# Logical scripts contain NAs error. 





# 