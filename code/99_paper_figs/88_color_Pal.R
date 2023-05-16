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
  theme(legend.position = "none") +
  labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")

TSNE_facet <- plotReducedDim(sce, dimred = "TSNE") +
  geom_point(aes(color = sce$bulkTypeSepHb)) +
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
    ), 
    size = 4,
    position = position_stack(vjust = 0.4),
    angle = -90,
    fontface = "bold",
    colour = "black",
    check_overlap = TRUE
  ) +
  theme_bw() +
  theme(legend.position = "None", 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(size = 16, hjust = 0.5)) +
  ggtitle("Bulk Deconvolution") +
  ylab("Proportion")

########## CHANGE PAL FUNCTION #################################################
change_Colors <- function(new_col_pal, title ){
  plot.top_1 <- plot_grid(TSNE + scale_colour_manual(values = new_col_pal), 
                       TSNE_facet + scale_colour_manual(values = new_col_pal), 
                       nrow = 1)
  
  plot.top_2 <- plot_grid(plot.top_1, 
                          comp_plot + scale_fill_manual(values = new_col_pal), 
                          nrow = 1,
                          rel_widths = c(.75, .25))
  
  plot_grid( plot.top_2,
             deconvo_plot + scale_fill_manual(values = alpha(new_col_pal, 0.8)),
             ncol = 1,
             labels = title
  ) 
}


### color pal test
# dr. mayanrd's ask
new_col_pal <- c("Oligo" = c("#4d5802"), 
                    "OPC"= c("#d3c871"), 
                    "OPC_noisy" = c("#A9A9A9"),
                    "Microglia" = c("#1c0f77"), 
                    "Astrocyte" = c("#8d363c"), 
                    "Endo" = c("#ee6c14"), 
                    "Excit.Neuron" = c("#71797E"), 
                    "Inhib.Thal" = c('#b5a2ff'),  
                    "Excit.Thal" = c("#9e4ad1"), 
                    "LHb" = c("#0085af"),
                    "MHb" = c("#fa246a")
)

# louise's preference 
new_col_pal_2 <- c("Oligo" = c("#d3c871"), 
                 "OPC"= c("#4d5802"), 
                 "OPC_noisy" = c("#A9A9A9"),
                 "Microglia" = c("#1c0f77"), 
                 "Astrocyte" = c("#8d363c"), 
                 "Endo" = c("#ee6c14"), 
                 "Excit.Neuron" = c("#71797E"), 
                 "Inhib.Thal" = c("#9e4ad1"),  
                 "Excit.Thal" = c('#b5a2ff'), 
                 "LHb" = c("#0085af"),
                 "MHb" = c("#fa246a")
)

# my idea : BULK COLOR PALATTE WINNER!!!!! #####################################
new_col_pal_3 <- c("Oligo" = c("#4d5802"), 
                 "OPC"= c("#d3c871"), 
                 "OPC_noisy" = c("#A9A9A9"),
                 "Microglia" = c("#222222"), 
                 "Astrocyte" = c("#8d363c"), 
                 "Endo" = c("#ee6c14"), 
                 "Excit.Neuron" = c("#71797E"), 
                 "Inhib.Thal" = c('#b5a2ff'),  
                 "Excit.Thal" = c("#9e4ad1"), 
                 "LHb" = c("#0085af"),
                 "MHb" = c("#fa246a")
)

################################################################################
pdf(file = here(plot_dir, "test.pdf"), width = 15, height = 12)
  change_Colors(new_col_pal, "Scheme 1")
dev.off()

pdf(file = here(plot_dir, "test_2.pdf"), width = 15, height = 12)
  change_Colors(new_col_pal_2, "Scheme 2")
dev.off()

pdf(file = here(plot_dir, "test_3.pdf"), width = 15, height = 12)
  change_Colors(new_col_pal_3, "Scheme 3")
dev.off()

#### Single Cell Resolution Color Trials #######################################
sn_test <- function(sn_color_pal, title) {
  plot1 <- plotReducedDim(sce, dimred = "TSNE") +
    geom_point(aes(color = sce$final_Annotations)) +
    scale_colour_manual(values = sn_color_pal) +
    theme(legend.position = "none") +
    labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")
  
  plot2 <- plotReducedDim(sce, dimred = "TSNE") +
    geom_point(aes(color = sce$final_Annotations)) +
    scale_colour_manual(values = sn_color_pal) +
    facet_wrap(~ sce$final_Annotations) +
    guides(color = guide_legend(title="Cell Type")) +
    labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")
  
  plot_grid(plot1, plot2,
            nrow = 1,
            labels = title)
} 

# Colors 1
sn_colors_1 <- c("Oligo" = c("#4d5802"), 
                   "OPC"= c("#d3c871"), 
                   "OPC_noisy" = c("#CCCCCC"),
                   "Microglia" = c("#222222"), 
                   "Astrocyte" = c("#8d363c"), 
                   "Endo" = c("#ee6c14"), 
                   "Excit.Neuron" = c("#666666"), 
                   "Inhib.Thal" = c('#b5a2ff'),  
                   "Excit.Thal" = c("#9e4ad1"),
                   "LHb.1" = c("#0085af"),
                   "LHb.2" = c("#0096FF"), 
                   "LHb.3" = c ("#89CFF0"), 
                   "LHb.4" = c("#6F8FAF"), 
                   "LHb.5" = c("#40E0D0"), 
                   "LHb.6" = c("#008080"),  
                   "LHb.7" = c("#7DF9FF"), 
                   "MHb.1" = c("#FF00FF"), 
                   "MHb.2" = c("#FAA0A0"),
                   "MHb.3" = c("#fa246a") 

)

pdf(file = here(plot_dir, "sn_test_1.pdf"), height = 7, width = 13)
  sn_test(sn_colors_1, "snRNA Scheme 1")
dev.off()


# Done. 