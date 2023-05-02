## May 2, 2023 - Bukola Ajanaku
# Using PCs, UMAPs, and TSNEs for supplement fig for thal vs habenula.
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
# Astrocyte       Endo Excit.Thal Inhib.Thal  LHb.1      LHb.2      LHb.3 
# 538         38       1800       7612        201        266        134 
# LHb.4      LHb.5      LHb.6      LHb.7      MHb.1      MHb.2      MHb.3 
# 477         83         39       1014        152        540         18 
# Microglia      Oligo        OPC 
# 145       2178       1796 

# creating plot directory
plot_dir <- here("plots", "99_paper_figs", "07_Hb_vs_Thal_Plots")
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

##### PLOTTING  ###########################################################
# post-harmony PC2 vs PC3
pdf(here(plot_dir, "PC2_vs_PC3_splitbyClusters.pdf"),
    width = 15, height = 8)
plot1 <- plotReducedDim(sce, dimred = "PCA", ncomponents = 2:3) +
  geom_point(aes(color = sce$final_Annotations)) +
  scale_colour_manual(values = cluster_colors) +
  theme(legend.position = "none")

plot2 <- plotReducedDim(sce, dimred = "PCA", ncomponents = 2:3) +
  geom_point(aes(color = sce$final_Annotations)) +
  scale_colour_manual(values = cluster_colors) +
  facet_wrap(~ sce$final_Annotations) +
  guides(color = guide_legend(title="Cell Type"))

plot_grid(plot1, plot2)
dev.off()


# Done.