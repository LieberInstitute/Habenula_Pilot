## April 11, 2023 - Bukola Ajanaku
# Exploring bisque results!
# qrsh -l mem_free=20G,h_vmem=20G

library(SummarizedExperiment)
library(here)
library(SingleCellExperiment)
library(DeconvoBuddies)
library(dplyr)

# loading sce object (sn data) post drop of Hb cluster! 
load(here("processed-data", "06_deconvolution", "sce_objects", "sce_first_bulkTypes.RDATA"),
     verbose = TRUE)
    # sce 

# loading cleaned rse object (bulk data)
load(here("processed-data", "02_bulk_qc", "count_data_bukola", 
          "rse_gene_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata"), verbose = TRUE)
    # rse_gene

# loading est_prop output from Bisque
load(here("processed-data", "06_deconvolution", "run_Bisque", 
     "est_prop_split_Hb_annotations.RDATA"), verbose = TRUE)

# grabbing relevant phenotype info for bulk data
pd <- colData(rse_gene) |>
  as.data.frame() |>
  select(Sample = RNum, Sex, PrimaryDx)
