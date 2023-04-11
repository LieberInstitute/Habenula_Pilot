## April 11, 2023 - Bukola Ajanaku
# Working on bulkRNA-seq deconvolution 
# qrsh -l mem_free=20G,h_vmem=20G

library(SummarizedExperiment)
library(here)
library(DeconvoBuddies)
library(SingleCellExperiment)

# loading sce object (sn data) post drop of Hb cluster! 
load(here("processed-data", "06_deconvolution", "sce_objects", "sce_first_bulkTypes.RDATA"),
     verbose = TRUE)
      # sce 
# loading cleaned rse object (bulk data)
load(here("processed-data", "02_bulk_qc", "count_data_bukola", 
          "rse_gene_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata"), verbose = TRUE)
      # rse_gene

###### CHECKING DATA ###########################################################
table(sce$bulkTypeSepHb)
    # Astrocyte   Endo Excit.Thal   Inhib.Thal   LHb        MHb  Microglia 
    # 538         38       1800       7612       2214        710        145 
    # Oligo       OPC 
    # 2178       1796


# Creating mean_ratios based on our specified annotations
ratios <- get_mean_ratio2(sce,
                          cellType_col = "cellType",
                          assay_name = "logcounts",
                          add_symbol = TRUE)
