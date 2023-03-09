## March 8, 2023 - Bukola Ajanaku
# Updating snAnno and replotting to verify identities
# qrsh -l mem_free=20G,h_vmem=20G

library("SingleCellExperiment")
library("here")
library("tidyverse")

# loading original sce object with just snAnno and other minor annotations
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
          "sce_post_09_clustered_qc.Rdata"))

# creating snAnno2 which will contain the combined MHb groups
sce$snAnno2 <- sce$snAnno

  ## LHb.6 is actually Endothelial. Total LHb is now 7 from 8.
sce$snAnno2[sce$snAnno2 == "LHB.6"] <- "Endo"

