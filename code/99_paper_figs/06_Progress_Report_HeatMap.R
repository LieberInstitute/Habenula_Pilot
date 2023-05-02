## May 2, 2023 - Bukola Ajanaku
# Remaking the Progress Report Heatmap with the cleaned OPC class! 
# qrsh -l mem_free=30G,h_vmem=30G

library("SingleCellExperiment")
library("here")
library("sessioninfo")
library("jaffelab")
library("ggplot2")
library("tidyverse")
library("ComplexHeatmap")
library("spatialLIBD")

# loading final sce object (no ambig cluster and can easily drop OPC noisy)
load(file = here("processed-data", "99_paper_figs", "sce_objects", 
                 "official_final_sce.RDATA"))

# creating plot_dir
plot_dir <- here("plots", "99_paper_figs", "06_Progress_Report_HeatMap")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

# dropping OPC_noisy
sce <- sce[ , which(sce$OPC_clean == "Yes")]

# Pseudobulking to create compressed sce object
## faking the pseudobulking function out by setting sample as all same sample
sce$FakeSample <- "Br1011"

set.seed(20220907) 
pb_sce <- registration_pseudobulk(sce, "final_Annotations", "FakeSample")

# list of marker genes 
official_markers = list(
  "Oligo" = c("MOBP"),
  "OPC" = c("PDGFRA"),
  "Micro" = c("CSF1R"),
  "Astro" = c("AQP4"),
  "Endo" = c("ITIH5"),
  "Thal" = c("LYPD6B"),
  "LHb.A" = c("LINC02653"), #  , ATP8B1
  "LHb.B" = c("AC073071.1"),
  "LHb.C" = c ("ENTHD1"),
  "LHb.D" = c("TLE2"),
  "LHb.E" = c("LINC01619"),
  "LHb.F" = c("TACR3"),
  "LHb.G" = c("AC008619.1"),
  "MHb.A" = c("EXOC1L"), 
  "MHb.B" = c("CHAT"),
  "MHb.C" = c("BHLHE22"),
  'Hb' = c("POU4F1"), # BARHL1
  "MHb" = c("CHRNB4"),
  "LHb" = c("HTR2C"),
  'Neuron' = c('SYT1'),
  'Exc_Neuron' = c('SLC17A6'), 
  'Inh_Neuron' = c('GAD1')
)

# cluster identities
row_namers <- c("Oligo",
                "OPC",
                "Microglia",
                "Astrocyte",
                "Endo",
                'Inhib.Thal', 
                'Excit.Thal', 
                "LHb.1",
                "LHb.2",
                "LHb.3",
                "LHb.4",
                "LHb.5",
                "LHb.6",
                "LHb.7",
                "MHb.1", 
                "MHb.2",
                "MHb.3"
)






# 