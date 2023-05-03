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
pd$Prop <- NA

# creating function that summarizes data
for(i in unique(sce$Sample)){
  pD <- pd[pd$Sample == i, ]
  
  # grab the total number of nuclei per Sample
  total_nuclei <- nrow(pD)
  
  # running through each cluster
  for(n in unique(pD$final_Annotations)){
    pD <- pD[pD$final_Annotations == n, ]
    
    # calculating proportion for given sample in respective cluster
    pD$Prop <- (nrow(pD) / total_nuclei) * 100
  }
  
  return(pD)
}

table(sce$final_Annotations)