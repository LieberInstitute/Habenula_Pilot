## April 17, 2023 - Bukola Ajanaku
# Plotting same TSNE plots but pre-harmony.
# qrsh -l mem_free=50G,h_vmem=50G

library("SingleCellExperiment")
library("here")
library("tidyverse")
library("scater")
library("ggplot2")
library("cowplot")

# loading old sce object (post qc sce object)
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
         "sce_uncorrected_PCA.Rdata"), verbose = TRUE)
# sce_uncorrected 

dim(sce_uncorrected)
# [1] 33848 17082

# loading official sce object 
load(here("processed-data", "99_paper_figs", "sce_objects", 
          "official_final_sce.RDATA"), verbose = TRUE)
# sce

dim(sce)
# [1] 33848 17031

# making sure colnames of sce_uncorrected are unique 
colnames(sce_uncorrected) <- paste0(sce_uncorrected$Sample, "_", sce_uncorrected$Barcode)

# subsetting sce_uncorrected to only the nuclei we've kept in sce
sce_uncorrected_clean <- sce_uncorrected[, which(colnames(sce_uncorrected) %in% colnames(sce))]

dim(sce_uncorrected_clean)
# [1] 33848 17031

# Now combining pheno data to sce_uncorrected_clean before saving.
colData(sce_uncorrected_clean) <- left_join(colData(sce_uncorrected_clean), 
                                  colData(sce))







