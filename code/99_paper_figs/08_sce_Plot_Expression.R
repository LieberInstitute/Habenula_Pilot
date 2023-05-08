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

#### get proportions before dropping ambig #####################################
colData(sce)[, c("final_Annotations", "Sample", "NeuN")]

prop_dirty <- as.data.frame(colData(sce)[, 
                                         c("final_Annotations", "Sample", "NeuN")]) |>
  group_by(Sample, final_Annotations, NeuN) |>
  summarize(n = n()) |>
  group_by(Sample) |>
  mutate(prop = n / sum(n))

# dropping the clusters we dropped
sce <- sce[, sce$final_Annotations != "OPC_noisy"]
sce <- sce[, sce$final_Annotations != "Excit.Neuron"]

#### proportions of nuclei using post-drop information #########################
pd <- as.data.frame(colData(sce))
table(pd$final_Annotations)
# Astrocyte   Endo Excit.Thal Inhib.Thal      LHb.1      LHb.2      LHb.3 
# 538         38       1800       7612        201        266        134 
# LHb.4      LHb.5      LHb.6      LHb.7      MHb.1      MHb.2      MHb.3 
# 477         83         39       1014        152        540         18 
# Microglia Oligo        OPC 
# 145       2178       1202 

prop_clean <- pd[,c("final_Annotations", "Sample", "NeuN")] |>
  group_by(Sample, final_Annotations, NeuN) |>
  summarize(n = n()) |>
  group_by(Sample) |>
  mutate(prop = n / sum(n))


### combines prop_dirty and prop_clean
prop_ambig_plus <- prop_dirty |>
  mutate(ambig = "Pre-drop") |>
  bind_rows(prop_clean |> mutate(ambig = "Post-drop"))

# plots composition plot using prop_clean and prop_dirty
comp_plot_both <- ggplot(data = prop_ambig_plus, aes(x = Sample, y = prop, fill = final_Annotations)) +
  geom_bar(stat = "identity") +
  geom_text(
    aes(
      label = ifelse(prop > 0.02, format(round(prop, 3), 3), ""),
    ),
    size = 2.5,
    position = position_stack(vjust = 0.5),
    color = "black"
  ) +
  labs(y = "Proportion", fill = "Cell Type") +
  facet_grid(fct_rev(ambig) ~ NeuN, scales = "free", space = "free") +
  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color = FALSE, fill = guide_legend(ncol = 1))

pdf(file = here(plot_dir, "full_Comp_Express_Plot.pdf"))
comp_plot_both 
dev.off()


#  