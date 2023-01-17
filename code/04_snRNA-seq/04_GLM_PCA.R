# January 17, 2023 - Bukola Ajanaku
# Creating pca plots for filtered (removed empty droplets), quality controlled 
# (dropped high mito, low library size, low detected features, and any genes 
# with 0 counts across samples) sce object.
# Based on:
# https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/code/03_build_sce/04_GLM_PCA.R
# qrsh -l mem_free = 50G, h_vmem = 50G

library("SingleCellExperiment")
library("scran")
library("scater")
library("scry")
library("here")
library("sessioninfo")

# loading post qc sce object 
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
           "sce_hb_postQC.Rdata"))
sce <- sce_hb_postQC
rm(sce_hb_postQC)

# Deviance featuring selection
set.seed(1234)

sce <- devianceFeatureSelection(sce,
        assay = "counts", fam = "binomial", sorted = F,
        batch = as.factor(sce$Sample))

 