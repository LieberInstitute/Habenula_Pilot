## May 10, 2023 - Bukola Ajanaku
# Making stepwise plot using cleaned TSNE, cleaned num nuclei plots, and cleaning
# composition plots.
# qrsh -l mem_free=50G,h_vmem=50G

# loading relevant libraries
library(SummarizedExperiment)
library(here)
library(SingleCellExperiment)
library(DeconvoBuddies)
library(tidyverse)
library(tibble)
library(cowplot)
library(scater)

# loading sce object with dropped ambig cluster
load(file = here("processed-data", "99_paper_figs", "sce_objects", 
                 "sce_final_preHbdrop.RDATA"), verbose = TRUE)
table(sce$final_Annotations)

# dropping Excit.Neuron and OPC_noisy clusters
sce <- sce[, sce$final_Annotations != "OPC_noisy"]
sce <- sce[, sce$final_Annotations != "Excit.Neuron"]

# creating plot_dir
plot_dir <- here("plots", "99_paper_figs", "09_step_Wise_Plots")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

# setting manual colors
cluster_colors <- c("Oligo" = c("#4d5802"), 
                         "OPC"= c("#9e4ad1"), 
                         "OPC_noisy" = c("#A9A9A9"),
                         "Microglia" = c("#1c0f77"), 
                         "Astrocyte" = c("#8d363c"), 
                         "Endo" = c("#ee6c14"), 
                         "Excit.Neuron" = c("#71797E"), 
                         "Inhib.Thal" = c("#d3c871"),  
                         "Excit.Thal" = c('#b5a2ff'), 
                         "LHb" = c("#0085af"),
                         "MHb" = c("#fa246a")
)

# grabbing bulk annotations 
sce$bulkTypeSepHb <- sce$final_Annotations

# Combining into 5 glia, two thalamus, 1 broad LHb, and 1 broad MHb.
sce$bulkTypeSepHb[sce$bulkTypeSepHb %in% grep("^LHb\\.", unique(sce$bulkTypeSepHb), value = TRUE)] <-  'LHb'
sce$bulkTypeSepHb[sce$bulkTypeSepHb %in% grep("^MHb\\.", unique(sce$bulkTypeSepHb), value = TRUE)] <-  'MHb'

# check levels
table(sce$bulkTypeSepHb)

# cleaned TSNE with facet_wrap
TSNE <- plotReducedDim(sce, dimred = "TSNE") +
  geom_point(aes(color = sce$bulkTypeSepHb)) +
  scale_colour_manual(values = cluster_colors) +
  theme(legend.position = "none")

TSNE_facet <- plotReducedDim(sce, dimred = "TSNE") +
  geom_point(aes(color = sce$bulkTypeSepHb)) +
  scale_colour_manual(values = cluster_colors) +
  facet_wrap(~ sce$bulkTypeSepHb) +
  guides(color = guide_legend(title="Cell Type"))

png(file = here(plot_dir, "bulk_clean_TSNE.png"), width = 900, height = 700)
  plot_grid(TSNE, TSNE_facet)
dev.off()

# number of nuclei per cell type post drop
num_nuc <- as.data.frame(colData(sce)[,c("final_Annotations", "bulkTypeSepHb", "Sample", "NeuN")]) |>
  group_by(Sample, bulkTypeSepHb, NeuN) |>
  mutate(n_nuc = n())

num_nuc_comp_plot <- ggplot(num_nuc, aes(x = final_Annotations, y = n_nuc, fill = bulkTypeSepHb)) +
  geom_col() +
  geom_text(aes(label = n_nuc), size = 2.5) +
  scale_fill_manual(values = cluster_colors) +
  theme_bw() +
  labs(y = "Number of Nuclei") 

png(file = here(plot_dir, "num_nuclei_post_clean.png"), width = 900, height = 700)
  num_nuc_comp_plot
dev.off()



# 