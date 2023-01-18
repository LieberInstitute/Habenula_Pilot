# January 18, 2023 - Bukola Ajanaku
# Normalizing sce object by different metrics and then plotting PCs, TSNEs, and
# UMAPs.
# Based on:
# 1) https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/code/03_build_sce/05_harmony_correction.R
# 2) https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/code/03_build_sce/05.5_harmony_correction_plots.R
# qrsh -l mem_free=50G,h_vmem=50G

# Loading relevant libraries
library("SingleCellExperiment")
library("harmony")
library("scater")
library("here")
library("sessioninfo")

# Loading sce_uncorrected
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
          "sce_uncorrected_PCA.Rdata"))

####### RUNNING HARMONY ########################################################
# Note: We used GLM PCA not normal PCA but harmony function searches for something
# named "PCA". Thus, we must submit GLM under the name PCA. After running harmony,
# we can do away with our fake PCA:

# Submitting GLM PCA under the name PCA:
reducedDim(sce_uncorrected, "PCA") <- reducedDim(sce_uncorrected, "GLMPCA_approx")

# Running harmony
sce_corrbySamp <- RunHarmony(sce_uncorrected, group.by.vars = "Sample", verbose = TRUE)
# sce_corrbyRun <- RunHarmony(sce_uncorrected, group.by.vars = "Run", verbose = TRUE)
    # Error in harmonyObj$init_cluster_cpp(0) : 
    # element-wise multiplication: incompatible matrix dimensions: 100x3 and 100x1

# Removing our redudant fake PCA name and keeping it under "GLMPCA_approx"
reducedDim(sce_uncorrected, "PCA") <- NULL

####### TSNE & UMAP ############################################################
set.seed(777)

## sce corrected by Sample
sce_corrbySamp <- runTSNE(sce_corrbySamp, dimred = "HARMONY")
sce_corrbySamp <- runUMAP(sce_corrbySamp, dimred = "HARMONY")

## sce corrected by Run
# sce_corrbyRun <- runTSNE(sce_corrbyRun, dimred = "HARMONY")
# sce_corrbyRun <- runUMAP(sce_corrbyRun, dimred = "HARMONY")

####### PLOTTING ###############################################################





####### saving ###############################################################






# 

