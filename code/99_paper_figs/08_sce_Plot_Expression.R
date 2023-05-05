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

# creating function that summarizes data
for(i in unique(sce$Sample)){
  pd[pd$Sample == i, ]$total_nuclei <- nrow(pd[pd$Sample == i, ])
  
  paste("Sample:", i)
  
  for(n in unique(pd$final_Annotations)){
    
    if(is.na(nrow(pd[pd$final_Annotations, ]))){
      tot_ct = 0
    } else{
      tot_ct = nrow(pd[pd$final_Annotations, ])
    }
    
    pd[pd$final_Annotations, ]$ct_nuclei <- tot_ct
    
    paste(n)
  }
}

# running through each cluster
for(n in unique(pd$final_Annotations)){
    pD <- pD[pD$ct_nuclei == n, ]
    ct_nuclei <- nrow(pD)
    
    # calculating proportion for given sample in respective cluster
    pD$Prop <- (ct_nuclei / total_nuclei) * 100
  }
}

table(sce$final_Annotations)