## April 11, 2023 - Bukola Ajanaku
# Exploring bisque results!
# qrsh -l mem_free=20G,h_vmem=20G

library(SummarizedExperiment)
library(here)
library(SingleCellExperiment)
library(DeconvoBuddies)
library(dplyr)
library(tibble)

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

# adding plot_dir
plot_dir <- here("plots", "06_deconvolution", "04_explore_Bisque")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

# grabbing relevant phenotype info for bulk data
pd <- colData(rse_gene) |>
  as.data.frame() |>
  select(Sample = RNum, BrNum, PrimaryDx)

prop_long <- est_prop$bulk.props |>
  as.data.frame() |>
  rownames_to_column("cellType") |>
  tidyr::pivot_longer(!cellType, names_to = "Sample", values_to = "prop") |>
  left_join(pd)
  
## create composition bar plots
pdf(here(plot_dir, "composition_bar_plot.pdf"), width = 20, height = 10)
  plot_composition_bar(prop_long = prop_long, sample_col = "Sample", x_col = "BrNum", ct_col = "cellType")
dev.off()
