## May 2, 2023 - Bukola Ajanaku
# Plotting cell-type expression pre and post drop per sample 
# qrsh -l mem_free=50G,h_vmem=50G

# loading relevant libraries
library(SummarizedExperiment)
library(here)
library(SingleCellExperiment)
library(DeconvoBuddies)
library(tidyverse)
library(tibble)

# loading sce object with dropped ambig cluster
load(file = here("processed-data", "99_paper_figs", "sce_objects", 
                 "sce_final_preHbdrop.RDATA"), verbose = TRUE)
table(sce$final_Annotations)
# Astrocyte         Endo Excit.Neuron   Excit.Thal   Inhib.Thal        LHb.1 
# 538           38           51         1800         7612          201 
# LHb.2        LHb.3        LHb.4        LHb.5        LHb.6        LHb.7 
# 266          134          477           83           39         1014 
# MHb.1        MHb.2        MHb.3    Microglia        Oligo          OPC 
# 152          540           18          145         2178         1202 
# OPC_noisy 
# 594 

# creating plot_dir
plot_dir <- here("plots", "99_paper_figs", "08_sce_Plot_Expression")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

# creating cell-type colors 
# adding color pallete (same color scheme used for progress report heatmap)
cluster_colors <- c( "Oligo" = c("#ba6425"), 
                     "OPC"= c("#987020"), 
                     "OPC_noisy" = c("#bf9146"),
                     "Microglia" = c("#5e0056"), 
                     "Astrocyte" = c("#e2693e"), 
                     "Endo" = c("#e295de"), 
                     "Excit.Neuron" = c("#709438"), 
                     "Inhib.Thal" = c("#eed7a1"),  
                     "Excit.Thal" = c('#E1F8DC'), 
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

####### FINAL ANNOTATIONS LEVEL #################################################
#### get proportions before dropping ambig #####################################
# grabbing proportion information
prop_dirty_sn <- as.data.frame(colData(sce)[, 
                                         c("final_Annotations", "Sample", "NeuN")]) |>
  group_by(Sample, final_Annotations, NeuN) |>
  summarize(n = n()) |>
  group_by(Sample) |>
  mutate(prop = n / sum(n))

# dropping the clusters we dropped
sce_drop <- sce[, sce$final_Annotations != "OPC_noisy"]
sce_drop <- sce_drop[, sce_drop$final_Annotations != "Excit.Neuron"]

#### proportions of nuclei using post-drop information #########################
pd_sn <- as.data.frame(colData(sce_drop))
table(pd_sn$final_Annotations)
# Astrocyte   Endo Excit.Thal Inhib.Thal      LHb.1      LHb.2      LHb.3 
# 538         38       1800       7612        201        266        134 
# LHb.4      LHb.5      LHb.6      LHb.7      MHb.1      MHb.2      MHb.3 
# 477         83         39       1014        152        540         18 
# Microglia Oligo        OPC 
# 145       2178       1202 

prop_clean_sn <- pd_sn[,c("final_Annotations", "Sample", "NeuN")] |>
  group_by(Sample, final_Annotations, NeuN) |>
  summarize(n = n()) |>
  group_by(Sample) |>
  mutate(prop = n / sum(n))

### combines prop_dirty and prop_clean
prop_ambig_plus_sn <- prop_dirty_sn |>
  mutate(ambig = "Pre-drop") |>
  bind_rows(prop_clean_sn |> mutate(ambig = "Post-drop"))

# plots composition plot using prop_clean and prop_dirty
comp_plot_both_sn <- ggplot(data = prop_ambig_plus_sn, aes(x = Sample, y = prop, fill = final_Annotations)) +
  geom_bar(stat = "identity") +
  geom_text(
    aes(
      label = ifelse(prop > 0.02, format(round(prop, 3), 3), "")
    ),
    size = 2.5,
    position = position_stack(vjust = 0.5),
    color = "black"
  ) +
  scale_fill_manual(values = c(cluster_colors)) +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
  labs(y = "Proportion", fill = "Cell Type") +
  facet_grid(ambig ~ NeuN, scales = "free", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color = FALSE, fill = guide_legend(ncol = 1))

pdf(file = here(plot_dir, "full_Comp_Express_Plot_finalAnnoLEVEL2.pdf"), width = 7, height = 11)
  comp_plot_both_sn
dev.off()

####### BULK COLLAPSE LEVEL ####################################################
cluster_colors_bulk <- c("Oligo" = c("#4d5802"), 
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

# creating bulk annotations level
sce$bulkTypeSepHb <- sce$final_Annotations

# Combining into 5 glia, two thalamus, 1 broad LHb, and 1 broad MHb.
sce$bulkTypeSepHb[sce$bulkTypeSepHb %in% grep("^LHb\\.", unique(sce$bulkTypeSepHb), value = TRUE)] <-  'LHb'
sce$bulkTypeSepHb[sce$bulkTypeSepHb %in% grep("^MHb\\.", unique(sce$bulkTypeSepHb), value = TRUE)] <-  'MHb'

# check levels
table(sce$bulkTypeSepHb)
    # Astrocyte         Endo Excit.Neuron   Excit.Thal   Inhib.Thal          LHb 
    # 538           38           51         1800         7612         2214 
    # MHb    Microglia        Oligo          OPC    OPC_noisy 
    # 710          145         2178         1202          594

#### get proportions before dropping ambig #####################################
prop_dirty_bulk <- as.data.frame(colData(sce)[, 
                                c("bulkTypeSepHb", "Sample", "NeuN")]) |>
  group_by(Sample, bulkTypeSepHb, NeuN) |>
  summarize(n = n()) |>
  group_by(Sample) |>
  mutate(prop = n / sum(n))


#### proportions of nuclei using post-drop information #########################
pd_bulk <- as.data.frame(colData(sce_drop))
table(pd_bulk$bulkTypeSepHb)
    # Astrocyte       Endo Excit.Thal Inhib.Thal        LHb        MHb  Microglia 
    # 538         38       1800       7612       2214        710        145 
    # Oligo        OPC 
    # 2178       1202

prop_clean_bulk <- pd_bulk[,c("bulkTypeSepHb", "Sample", "NeuN")] |>
  group_by(Sample, bulkTypeSepHb, NeuN) |>
  summarize(n = n()) |>
  group_by(Sample) |>
  mutate(prop = n / sum(n))

### combines prop_dirty and prop_clean
prop_ambig_plus_bulk <- prop_dirty_bulk |>
  mutate(ambig = "Pre-drop") |>
  bind_rows(prop_clean_bulk |> mutate(ambig = "Post-drop")) |> 
  mutate(ct_levels = factor(bulkTypeSepHb, levels = 
                              c("Excit.Neuron",
                                "Astrocyte", 
                                "Endo", 
                                "Microglia", 
                                "Oligo", 
                                "OPC_noisy",
                                "OPC",
                                "Inhib.Thal", 
                                "Excit.Thal" , 
                                "MHb", 
                                "LHb")) ) |>
  arrange(ct_levels)

# plots composition plot using prop_clean and prop_dirty
comp_plot_both_bulk <- ggplot(data = prop_ambig_plus_bulk, aes(x = Sample, 
                              y = prop, fill = bulkTypeSepHb)) +
  geom_bar(stat = "identity") +
  geom_text(
    aes(
      label = ifelse(prop > 0.02, format(round(prop, 3), 3), "")
    ),
    size = 2.5,
    position = position_stack(vjust = 0.5),
    color = "black"
  ) +
  scale_fill_manual(values = c(cluster_colors_bulk)) +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
  labs(y = "Proportion", fill = "Cell Type") +
  facet_grid(fct_rev(ambig) ~ NeuN, scales = "free", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color = FALSE, fill = guide_legend(ncol = 1))

pdf(file = here(plot_dir, "full_Comp_Express_Plot_bulkAnnoLEVEL_OFFICIAL.pdf"), width = 7, height = 11)
  comp_plot_both_bulk
dev.off()

# plotting total nuclei information per sample
barplot_n_nuc_bulk <- ggplot(prop_ambig_plus_bulk, 
  aes(x = Sample, y = n, fill = bulkTypeSepHb)) +
  geom_col() +
  geom_text(aes(label = n), size = 2.5) +
  scale_fill_manual(values = cluster_colors_bulk) +
  theme_bw() +
#  theme(legend.position = "None", axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
  labs(y = "Number of Nuclei") +
  facet_grid(fct_rev(ambig) ~ NeuN, scales = "free", space = "free")

pdf(file = here(plot_dir, "num_Nuc_Comp_Plot_bulkAnnoLEVEL.pdf"))
  barplot_n_nuc_bulk
dev.off()

# plotting total nuclei information per sample
sum_nuc_ambig_plus_prop <- prop_ambig_plus_bulk |>
  group_by(ambig, Sample, bulkTypeSepHb) |>
  summarize(n_across_samps = sum(n))

barplot_n_nuc_bulk_tot <- ggplot(sum_nuc_ambig_plus_prop, 
         aes(x = bulkTypeSepHb, y = n_across_samps, fill = bulkTypeSepHb)) +
  geom_col() +
  geom_text(aes(label = n_across_samps), size = 2.5) +
  scale_fill_manual(values = cluster_colors_bulk) +
  theme_bw() +
  #  theme(legend.position = "None", axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
  labs(y = "Number of Nuclei") +
  facet_wrap( ~ ambig, ncol = 1)

pdf(file = here(plot_dir, "num_Nuc_Comp_Plot_bulkAnnoLEVEL_overall.pdf"), width = 10, height = 9)
 barplot_n_nuc_bulk_tot
dev.off()



# 