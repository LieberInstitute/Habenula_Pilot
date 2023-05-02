## April 18, 2023 - Bukola Ajanaku
# Plotting same UMAP plots but pre-harmony.
# qrsh -l mem_free=50G,h_vmem=50G

library("SingleCellExperiment")
library("here")
library("tidyverse")
library("scater")
library("ggplot2")
library("cowplot")

# loading old sce object (post qc sce object)
load(here("processed-data", "99_paper_figs", "sce_objects", 
    "official_final_uncorrected_sce.RDATA"), verbose = TRUE)
# sce_uncorrected_clean

#### Prepping for plotting #####################################################
# creating plot directory
plot_dir <- here("plots", "99_paper_figs", "02_UMAPs", "Pre-Harmony")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

# adding color pallete (same color scheme used for progress report heatmap)
cluster_colors <- c( "Oligo" = c("#475c6c"), 
                     "OPC"= c("#899499"), 
                     "Microglia" = c("#bfb5b2"), 
                     "Astrocyte" = c("#c7bbc9"), 
                     "Endo" = c("#8a8583"), 
                     "Excit.Neuron" = c("#cd8b62"), 
                     "Inhib.Thal" = c("#eed7a1"),  
                     "Excit.Thal" = c('#f7efd2'), 
                     "LHb.1" = c("#00FFFF"),
                     "LHb.2" = c("#0096FF"), 
                     "LHb.3" = c ("#1434A4"), 
                     "LHb.4" = c("#00008B"), 
                     "LHb.5" = c("#40E0D0"), 
                     "LHb.6" = c("#008080"),  
                     "LHb.7" = c("#7DF9FF"), 
                     "MHb.1" = c("#800020"), 
                     "MHb.2" = c("#D70040"),
                     "MHb.3" = c("#D2042D") 
)

## create sce_scorted and unsorted based on NeuN
sce_unc_sorted <- sce_uncorrected_clean[, which(sce_uncorrected_clean$NeuN == "NeuN.Sorted")]
sce_unc_unsorted <- sce_uncorrected_clean[, which(sce_uncorrected_clean$NeuN == "NeuN.Unsorted")]

##### PLOTTING UMAPs ###########################################################
# Pre-Harmony UMAP colored and split by Sample
pdf(here(plot_dir, "UMAP_uncorrected_by_Sample.pdf"), width = 11, height = 7)
plot1 <-plotReducedDim(sce_uncorrected_clean, dimred = "UMAP") +
    geom_point(aes(color = sce_uncorrected_clean$Sample), alpha = 0.4) + 
    theme(legend.position = "none")

plot2 <-plotReducedDim(sce_uncorrected_clean, dimred = "UMAP") +
    geom_point(aes(color = sce_uncorrected_clean$Sample), alpha = 0.4) +
    facet_wrap(~ sce_uncorrected_clean$Sample) +
    guides(color = guide_legend(title="Sample ID"))

plot_grid(plot1, plot2)
dev.off()

# Pre-Harmony UMAP colored by Sample and split by Run
pdf(here(plot_dir, "UMAP_uncorrected_by_Sample_splitby_Run.pdf"), width = 11, height = 5)
plot1 <-plotReducedDim(sce_uncorrected_clean, dimred = "UMAP") +
  geom_point(aes(color = sce_uncorrected_clean$Sample), alpha = 0.4) + 
  theme(legend.position = "none")

plot2 <-plotReducedDim(sce_uncorrected_clean, dimred = "UMAP") +
  geom_point(aes(color = sce_uncorrected_clean$Sample), alpha = 0.4) +
  facet_wrap(~ sce_uncorrected_clean$Run) +
  guides(color = guide_legend(title="Sample ID"))

plot_grid(plot1, plot2, rel_widths = c(2,5))
dev.off()


# Done.

