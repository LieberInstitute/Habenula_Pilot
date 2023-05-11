## May 11, 2023 - Bukola Ajanaku
# Working on hockey stick plots and top 50 marker gene table using cleaned sce
# object (no OPC_noisy or Exit.Neuron) and at 9 cluster (bulk collapse) level.
# qrsh -l mem_free=50G,h_vmem=50G

library(SummarizedExperiment)
library(here)
library(DeconvoBuddies)
library(SingleCellExperiment)
library(jaffelab)
library(dplyr)
library(ggplot2)

# loading official sce object
load(file = here("processed-data", "99_paper_figs", "sce_objects", 
                 "official_final_sce.RDATA"), verbose = TRUE)

# has dropped Excit.Neuron and OPC_noisy
table(sce$final_Annotations)
  # Astrocyte       Endo Excit.Thal Inhib.Thal      LHb.1      LHb.2      LHb.3 
  # 538         38       1800       7612        201        266        134 
  # LHb.4      LHb.5      LHb.6      LHb.7      MHb.1      MHb.2      MHb.3 
  # 477         83         39       1014        152        540         18 
  # Microglia      Oligo        OPC 
  # 145       2178       1202 

###### Adding bulk collapsed annotations to sce object #########################
sce$bulkTypeSepHb <- sce$final_Annotations
# making separated Hb (2)
sce$bulkTypeSepHb[sce$bulkTypeSepHb %in% grep("^LHb\\.", unique(sce$bulkTypeSepHb), value = TRUE)] <-  'LHb'
sce$bulkTypeSepHb[sce$bulkTypeSepHb %in% grep("^MHb\\.", unique(sce$bulkTypeSepHb), value = TRUE)] <-  'MHb'

table(sce$bulkTypeSepHb)
  # Astrocyte       Endo Excit.Thal Inhib.Thal        LHb        MHb  Microglia 
  # 538         38       1800       7612       2214        710        145 
  # Oligo        OPC 
  # 2178       1202 











# 