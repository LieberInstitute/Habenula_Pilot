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
# Post-Harmonnization
pdf(here(plot_dir, "UMAP_harmony_by_finalAnno_PRE-HARMONY.pdf"))
plotReducedDim(sce_uncorrected_clean, dimred = "UMAP") +
  geom_point(aes(color = sce_uncorrected_clean$final_Annotations)) + 
  scale_colour_manual(values = cluster_colors) + 
  guides(color = guide_legend(title="Cell Type"))
dev.off()

pdf(here(plot_dir, "UMAP_harmony_by_final_Annotations_splitbyfinalAnno_PRE-HARMONY.pdf"), 
    width = 10)
plot1 <- plotReducedDim(sce_uncorrected_clean, dimred = "UMAP") +
  geom_point(aes(color = sce_uncorrected_clean$final_Annotations)) + 
  scale_colour_manual(values = cluster_colors) + 
  theme(legend.position = "none")

plot2 <- plotReducedDim(sce_uncorrected_clean, dimred = "UMAP") +
  geom_point(aes(color = sce_uncorrected_clean$final_Annotations)) + 
  scale_colour_manual(values = cluster_colors) +
  facet_wrap(~ sce_uncorrected_clean$final_Annotations) + 
  guides(color = guide_legend(title="Cell Type"))

plot_grid(plot1, plot2)
dev.off()

# for NeuN sorting and Non-NeuN sorting
plot_sorted <- plotReducedDim(sce_unc_sorted, dimred = "UMAP") +
  geom_point(aes(color = sce_unc_sorted$final_Annotations), alpha = 0.4) + 
  scale_colour_manual(values = cluster_colors) +
  facet_grid(sce_unc_sorted$NeuN ~ sce_unc_sorted$Sample) + 
  guides(color = guide_legend(title="Cell Type"))

plot_unsorted <-   plotReducedDim(sce_unc_unsorted, dimred = "UMAP") +
  geom_point(aes(color = sce_unc_unsorted$final_Annotations), alpha = 0.4) + 
  scale_colour_manual(values = cluster_colors) +
  facet_grid(sce_unc_unsorted$NeuN ~ sce_unc_unsorted$Sample) + 
  guides(color = guide_legend(title="Cell Type"))


pdf(here(plot_dir, "UMAP_harmony_by_finalAnno_splitbySampleAndSorting_PRE-HARMONY.pdf"), width = 13, height = 9)
plot_grid(
  plot_sorted,
  plot_unsorted,
  ncol = 1
)
dev.off()

pdf(here(plot_dir, "UMAP_harmony_by_final_Annotations_splitbyRun_PRE-HARMONY.pdf"), 
    width = 20, height = 8)
plot1 <- plotReducedDim(sce_uncorrected_clean, dimred = "UMAP") +
  geom_point(aes(color = sce_uncorrected_clean$final_Annotations)) + 
  scale_colour_manual(values = cluster_colors) +
  facet_wrap(~ sce_uncorrected_clean$Run) + 
  theme(legend.position = "none")

plot2 <- plotReducedDim(sce_uncorrected_clean, dimred = "UMAP") +
  geom_point(aes(color = sce_uncorrected_clean$final_Annotations)) + 
  scale_colour_manual(values = cluster_colors) +
  facet_wrap(~ sce_uncorrected_clean$final_Annotations) + 
  guides(color = guide_legend(title="Cell Type"))

plot_grid(plot1, plot2)
dev.off()


# Done.

