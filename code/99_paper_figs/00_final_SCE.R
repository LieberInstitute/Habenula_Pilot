## April 6, 2023 - Bukola Ajanaku
# TSNEs by my celltypes and custom color pallete 
# qrsh -l mem_free=50G,h_vmem=50G

library("SingleCellExperiment")
library("here")
library("tidyverse")
library("scater")
library("ggplot2")
library("cowplot")


# Loading sce object
load(here("processed-data", "99_paper_figs",  "sce_objects", "sce_final_preHbdrop.RDATA"),
     verbose = TRUE)
# sce_final_preHbdrop, sce_sorted, sce_unsorted

sce <- sce_final_preHbdrop
rm(sce_final_preHbdrop)

table(sce$final_Annotations)

# dropping redundant colData 
dropColData <- c("ct_Erik", "louvain10", "louvain20", "louvain50",
                 "wT_20_Erik", "wT_50_Erik", "groupErik", "wt10_ANNO",
                 "wt50_ANNO", "cellType_wT50", "cellType_wT10", 
                 "cellType_wT20", "snAnno", "splitProbClusts", 
                 "splitSNType", "snAnno2", "splitSNType2")

for(i in dropColData){
  colData(sce)[,i] <- NULL
}

# dropping small Hb cluster aka the Excit.Neuron Cluster
sce <- sce[, which(sce$final_Annotations != "Excit.Neuron")]

# adding OPC_clean column to sort for noisy OPC samples vs clean OPC samples 
# within the same cluster

# creating OPC_clean default
sce$OPC_clean <- "Yes"

# adding rownames of colData as a row for easier subsetting
sce$Rows <- rownames(colData(sce))

# OPC_noisy Samples = c("Br5555", "Br1204", "Br1092")
OPC_noisy = c("Br5555", "Br1204", "Br1092")

# grabbing barcodes for noisy OPC
RowNos <- sce[, sce[, which(sce$Sample %in% 
                  OPC_noisy)]$final_Annotations == "OPC"]$Rows

# adding Nos to RowNos 
sce[, sce$Rows %in% RowNos]$OPC_clean <- "No"

# check
table(sce$Sample, sce$OPC_clean)




# 

