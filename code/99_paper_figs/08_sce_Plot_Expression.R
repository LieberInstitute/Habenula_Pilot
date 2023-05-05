## May 2, 2023 - Bukola Ajanaku
# Plotting cell-type expression pre and post drop per sample 
# qrsh -l mem_free=20G,h_vmem=20G

# loading relevant libraries
library(SummarizedExperiment)
library(here)
library(SingleCellExperiment)
library(DeconvoBuddies)
library(tidyverse)

# loading sce object with dropped ambig cluster
load(file = here("processed-data", "99_paper_figs", "sce_objects", 
                 "official_final_sce.RDATA"))

# dropping OPC_noisy
table(sce$OPC_clean)
    # No   Yes 
    # 594 16437 
sce <- sce[, which(sce$OPC_clean == "Yes")]
# check again
table(sce$OPC_clean)
    # Yes 
    # 16437 

# grabbing relevant columns
pd <- colData(sce)[,c("Sample", "final_Annotations", "NeuN", "Run")]
pd$total_nuclei <- NA
pd$ct_nuclei <- NA

# creating function that collects nuclei data
for(i in unique(pd$Sample)){

  pd[pd$Sample == i, ]$total_nuclei <- nrow(pd[pd$Sample == i, ])

  for(p in unique(pd$final_Annotations)){
    
    if(nrow(pd[pd$Sample == i & pd$final_Annotations == p, ]) == 0){
      tot_ct = 0
    } else{
      tot_ct = nrow(pd[pd$Sample == i & pd$final_Annotations == p, ])
      pd[pd$Sample == i & pd$final_Annotations == p, ]$ct_nuclei <- tot_ct
    }
  }
}

pd$prop <- (pd$ct_nuclei / pd$total_nuclei)*100
# testing to make sure it adds up to %100
  # test <- pd[pd$Sample == i,]
  # tester <- unique(test$prop)
  # sum(tester)
   # [1] 99.90193








# 