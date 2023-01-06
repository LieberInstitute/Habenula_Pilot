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

# Loading raw data:
load(here("processed-data","09_snRNA-seq_re-processed","20220601_human_hb_processing.rda"))

# Grabbing locations of raw data for reading in
addressRawData <- unique(sce.all.hb$Sample)

# Reading in data
testRead <- read10xCounts(addressRawData[1], col.names=TRUE) 
