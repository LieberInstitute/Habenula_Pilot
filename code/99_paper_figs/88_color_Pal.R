## May 15, 2023 - Bukola Ajanaku
# Testing how different color palletes will look on our presenting data. 
# qrsh -l mem_free=20G,h_vmem=20G

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

# loading bulk deconvolution data 
load(file = here("processed-data", "99_paper_figs", "sce_objects", "prop_long.RDATA"),
     verbose = TRUE)

# creating plot_dir
plot_dir <- here("plots", "99_paper_figs", "88_color_Pal")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

# creating plot examples

######## TSNEs #################################################################
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
  theme(legend.position = "none") +
  labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")

TSNE_facet <- plotReducedDim(sce, dimred = "TSNE") +
  geom_point(aes(color = sce$bulkTypeSepHb)) +
  scale_colour_manual(values = cluster_colors) +
  facet_wrap(~ sce$bulkTypeSepHb) +
  guides(color = guide_legend(title="Cell Type"))

######### COMPOSITION PLOTs ####################################################
prop_df <- as.data.frame(colData(sce)[, c("bulkTypeSepHb", "Sample")]) |>
  group_by(Sample, bulkTypeSepHb) |>
  summarize(n = n()) |>
  group_by(Sample) |>
  mutate(prop = n / sum(n))

comp_plot <- ggplot(prop_df, 
                    aes(x = Sample, y = prop, fill = bulkTypeSepHb)) +
  geom_col() +
  geom_text(
    aes(
      label = ifelse(prop > 0.02, format(round(prop, 3), 3), "")
    ), 
    size = 3,
    position = position_stack(vjust = 0.5),
    color = "white",
  ) +
  scale_fill_manual(values = cluster_colors) +
  theme_bw() +
  theme(legend.position = "None", 
        axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
  labs(y = "Proportion")

########## BULK DECONVOLUTION PLOTs ############################################
deconvo_plot <- ggplot(prop_long, 
                    aes(x = Br_Order, y = prop_perc, fill = factor_CT)) +
  geom_col() +
  geom_text(
    data = subset(prop_long, prop_perc >= 1),
    aes(
      label = round(prop_perc, 1),
      check_overlap = TRUE
    ), 
    size = 4,
    position = position_stack(vjust = 0.4),
    angle = -90,
    fontface = "bold",
    colour = "black"
  ) +
  scale_fill_manual(values = alpha(color_bulk_clusters, 0.8)) +
  theme_bw() +
  theme(legend.position = "None", 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(size = 16, hjust = 0.5)) +
  ggtitle("Bulk Deconvolution") +
  ylab("Proportion")

########## CHANGE PAL FUNCTION #################################################

change_Colors <- function(new_col_pal){
  plot_grid( 
    TSNE + scale_colour_manual(values = new_col_pal), 
    TSNE_facet + scale_colour_manual(values = new_col_pal), 
    comp_plot + scale_fill_manual(values = cluster_colors),
    deconvo_plot + scale_fill_manual(values = alpha(color_bulk_clusters, 0.8))
  )
}



###





# 