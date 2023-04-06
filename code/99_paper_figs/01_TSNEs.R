## April 6, 2023 - Bukola Ajanaku
# TSNEs by my celltypes and custom color pallete 
# qrsh -l mem_free=20G,h_vmem=20G

library("SingleCellExperiment")
library("here")
library("tidyverse")


# Loading sce object with snAnno2 (Has combined MHb clusters and adjusted LHb clusters. 
# We probs wanna rename the general "Hb" cluster to something like Exc.Neuron!!)
load(here("processed-data", "04_snRNA-seq", "sce_objects", "sce_with_snAnno2.RDATA"), 
     verbose = TRUE)

# summary()