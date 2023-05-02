## April 6, 2023 - Bukola Ajanaku
# UMAPs by my celltypes and custom color pallete 
# qrsh -l mem_free=20G,h_vmem=20G

library("SingleCellExperiment")
library("here")
library("tidyverse")
library("scater")
library("ggplot2")
library("cowplot")


# Loading sce object
load(here("processed-data", "99_paper_figs", "sce_objects", 
          "official_final_sce.RDATA"), verbose = TRUE)
  # sce

table(sce$final_Annotations)
  # Astrocyte         Endo Excit.Neuron   Excit.Thal   Inhib.Thal        LHb.1 
  # 538           38           51         1800         7612          201 
  # LHb.2        LHb.3        LHb.4        LHb.5        LHb.6        LHb.7 
  # 266          134          477           83           39         1014 
  # MHb.1        MHb.2        MHb.3    Microglia        Oligo          OPC 
  # 152          540           18          145         2178         1796 

# creating plot directory
plot_dir <- here("plots", "99_paper_figs", "02_UMAPs", "Post-Harmony")
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
sce_sorted <- sce[, which(sce$NeuN == "NeuN.Sorted")]
sce_unsorted <- sce[, which(sce$NeuN == "NeuN.Unsorted")]


##### PLOTTING UMAPs ###########################################################
# Post-Harmony UMAP colored and split by Sample
pdf(here(plot_dir, "UMAP_harmony_by_Sample.pdf"), width = 11, height = 7)
plot1 <-plotReducedDim(sce, dimred = "UMAP") +
  geom_point(aes(color = sce$Sample), alpha = 0.4) + 
  theme(legend.position = "none")

plot2 <-plotReducedDim(sce, dimred = "UMAP") +
  geom_point(aes(color = sce$Sample), alpha = 0.4) +
  facet_wrap(~ sce$Sample) +
  guides(color = guide_legend(title="Sample ID"))

plot_grid(plot1, plot2)
dev.off()

# Post-Harmony UMAP colored by Sample and split by Run
pdf(here(plot_dir, "UMAP_harmony_by_Sample_splitby_Run.pdf"), width = 11, height = 5)
plot1 <-plotReducedDim(sce, dimred = "UMAP") +
  geom_point(aes(color = sce$Sample), alpha = 0.4) + 
  theme(legend.position = "none")

plot2 <-plotReducedDim(sce, dimred = "UMAP") +
  geom_point(aes(color = sce$Sample), alpha = 0.4) +
  facet_wrap(~ sce$Run) +
  guides(color = guide_legend(title="Sample ID"))

plot_grid(plot1, plot2, rel_widths = c(2,5))
dev.off()

# Post-Harmony UMAP colored by cell-type and split by cell-type
pdf(here(plot_dir, "UMAP_harmony_by_final_Annotations.pdf"), width = 11, height = 5)
plot1 <-plotReducedDim(sce, dimred = "UMAP") +
  geom_point(aes(color = sce$final_Annotations), alpha = 0.4) + 
  theme(legend.position = "none") + 
  scale_colour_manual(values = cluster_colors) 

plot2 <-plotReducedDim(sce, dimred = "UMAP") +
  geom_point(aes(color = sce$final_Annotations), alpha = 0.4) +
  facet_wrap(~ sce$final_Annotations) +
  guides(color = guide_legend(title="Sample ID")) + 
  scale_colour_manual(values = cluster_colors)

plot_grid(plot1, plot2)
dev.off()

# Done.