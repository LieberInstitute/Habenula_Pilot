## May 2, 2023 - Bukola Ajanaku
# Post QC sce object being shown per annotated cluster.
# qrsh -l mem_free=30G,h_vmem=30G

library("SingleCellExperiment")
library("here")
library("ggplot2")
library("sessioninfo")
library("scater")
library("cowplot")

# loading sce object post-cleaning and annotations (includes Hb cluster)
load(here("processed-data", "99_paper_figs",  "sce_objects", "sce_final_preHbdrop.RDATA"),
     verbose = TRUE)
# not needed 
sce_sorted <- NULL
sce_unsorted <- NULL

# renaming for ease
sce <- sce_final_preHbdrop
rm(sce_final_preHbdrop)

# sourcing official color palette 
source(file = here("code", "99_paper_figs", "source_colors.R"))
# bulk_colors and sn_colors

# creating plot_dir
plot_dir <- here("plots", "99_paper_figs", "05_QC_Boxplots")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

# plotting per metric on final_Annotations 
# also has OPC_noisy class
pdf(file = here(plot_dir, "sce_Annotated_with_Ambig_QC_boxplot.pdf"), width = 12, height = 12)

plot1 <- ggcells(sce, mapping = aes(x = final_Annotations, y = doubletScore, fill = final_Annotations)) +
      geom_boxplot() +
      scale_fill_manual(values = sn_colors) +
      theme_linedraw() +
      theme(
        legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) + labs (y = "Doublet Score")

plot2 <- ggcells(sce, mapping = aes(x = final_Annotations, y = subsets_Mito_percent, fill = final_Annotations)) +
  geom_boxplot() +
  scale_fill_manual(values = sn_colors) +
  theme_linedraw() +
  theme(
    legend.position = "none", axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) + labs (y = "Percent Mito")

plot3 <- ggcells(sce, mapping = aes(x = final_Annotations, y = sum, fill = final_Annotations)) +
  geom_boxplot() +
  scale_fill_manual(values = sn_colors) +
  theme_linedraw() +
  theme(
    legend.position = "none", axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) + labs (y = "Library Size")

plot4 <- ggcells(sce, mapping = aes(x = final_Annotations, y = detected, fill = final_Annotations)) +
  geom_boxplot() +
  scale_fill_manual(values = sn_colors) +
  theme_linedraw() +
  theme(
    legend.position = "none", axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) + labs (y = "Detected Features")

plot_grid(plot1, plot2, plot3, plot4,
          labels = c("A", "B", "C", "D"),
          ncol = 1)

dev.off()

# dropping Excit.Neuron
sce_drop <- sce[, which(!sce$final_Annotations == "Excit.Neuron")]
sce_drop <- sce[, which(!sce$final_Annotations == "OPC_noisy")]

# check 
table(sce_drop$final_Annotations)
# Astrocyte       Endo Excit.Thal Inhib.Thal      LHb.1      LHb.2      LHb.3 
# 538         38       1800       7612        201        266        134 
# LHb.4      LHb.5      LHb.6      LHb.7      MHb.1      MHb.2      MHb.3 
# 477         83         39       1014        152        540         18 
# Microglia      Oligo        OPC 
# 145       2178       1796

# plotting per metric on final_Annotations without ambig class
# also has OPC_noisy class
pdf(file = here(plot_dir, "sce_No_Ambig_Annotated_QC_boxplot.pdf"), width = 12, height = 12)

plot1 <- ggcells(sce_drop, mapping = aes(x = final_Annotations, y = doubletScore, fill = final_Annotations)) +
  geom_boxplot() +
  scale_fill_manual(values = sn_colors) +
  theme_linedraw() +
  theme(
    legend.position = "none", axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) + labs (y = "Doublet Score")

plot2 <- ggcells(sce_drop, mapping = aes(x = final_Annotations, y = subsets_Mito_percent, fill = final_Annotations)) +
  geom_boxplot() +
  scale_fill_manual(values = sn_colors) +
  theme_linedraw() +
  theme(
    legend.position = "none", axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) + labs (y = "Percent Mito")

plot3 <- ggcells(sce_drop, mapping = aes(x = final_Annotations, y = sum, fill = final_Annotations)) +
  geom_boxplot() +
  scale_fill_manual(values = sn_colors) +
  theme_linedraw() +
  theme(
    legend.position = "none", axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) + labs (y = "Library Size")

plot4 <- ggcells(sce_drop, mapping = aes(x = final_Annotations, y = detected, fill = final_Annotations)) +
  geom_boxplot() +
  scale_fill_manual(values = sn_colors) +
  theme_linedraw() +
  theme(
    legend.position = "none", axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) + labs (y = "Detected Features")

plot_grid(plot1, plot2, plot3, plot4,
          labels = c("A", "B", "C", "D"),
          ncol = 1)

dev.off()


# Done.