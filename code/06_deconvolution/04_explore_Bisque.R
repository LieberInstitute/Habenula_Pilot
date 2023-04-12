## April 11, 2023 - Bukola Ajanaku
# Exploring bisque results!
# qrsh -l mem_free=20G,h_vmem=20G

library(SummarizedExperiment)
library(here)
library(SingleCellExperiment)
library(DeconvoBuddies)
library(dplyr)
library(tibble)
library(ggplot2)

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

# custom color scheme
color_bulk_clusters <- 
  c( "Oligo" = c("#A9A9A9"), # dark grey
     "OPC"= c("#7393B3"), # blue grey
     "Microglia" = c("#E5E4E2"), # platinum
     "Astrocyte" = c("#36454F"), # ash grey
     "Endo" = c("#848884"), # smoke
     "Inhib.Thal" = c('#2AAA8A'), # jungle green
     "Excit.Thal" = c("#478778"), # lincoln green
     "LHb" = c("#DE3163"), # cerise
     "MHb" = c("#00FFFF") # aqua
)

# grabbing relevant phenotype info for bulk data
pd <- colData(rse_gene) |>
  as.data.frame() |>
  select(Sample = RNum, BrNum, PrimaryDx)

prop_long <- est_prop$bulk.props |>
  as.data.frame() |>
  rownames_to_column("cellType") |>
  tidyr::pivot_longer(!cellType, names_to = "Sample", values_to = "prop") |>
  left_join(pd) 

prop_long$HbSum <- NA

# creating a column to organize plot by Hb sum.
for(i in unique(prop_long$BrNum)) {
  prop_sub <- prop_long[prop_long$BrNum == i,] 
  HbSummer <- sum(prop_sub[prop_sub$cellType %in% c("MHb", "LHb"), ]$prop)
  prop_long[prop_long$BrNum == i,]$HbSum <- HbSummer
}

prop_long$Hb_Order <- order(prop_long$HbSum)

prop_long <- prop_long %>% arrange(desc(Hb_Order))

## create composition bar plots
pdf(here(plot_dir, "composition_bar_plot_bulkTypeSepHb_2.pdf"), width = 20, height = 10)
  plot_composition_bar(prop_long = prop_long, sample_col = "BrNum",
                        x_col = "BrNum", ct_col = "cellType") + 
    scale_fill_manual(values = color_bulk_clusters) +
    aes(x = reorder(prop_long$BrNum, prop_long$Hb_Order), y = prop_long$prop)
dev.off()






#
