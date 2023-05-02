## April 17, 2023 - Bukola Ajanaku
# Plotting same TSNE plots but pre-harmony.
# qrsh -l mem_free=50G,h_vmem=50G

library("SingleCellExperiment")
library("here")
library("tidyverse")
library("scater")
library("ggplot2")
library("cowplot")

# loading old sce object (post qc sce object)
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
         "sce_uncorrected_PCA.Rdata"), verbose = TRUE)
# sce_uncorrected 

dim(sce_uncorrected)
# [1] 33848 17082

# loading official sce object 
load(here("processed-data", "99_paper_figs", "sce_objects", 
          "official_final_sce.RDATA"), verbose = TRUE)
# sce

dim(sce)
# [1] 33848 17031

# making sure colnames of sce_uncorrected are unique 
colnames(sce_uncorrected) <- paste0(sce_uncorrected$Sample, "_", sce_uncorrected$Barcode)

# subsetting sce_uncorrected to only the nuclei we've kept in sce
sce_uncorrected_clean <- sce_uncorrected[, which(colnames(sce_uncorrected) %in% colnames(sce))]

dim(sce_uncorrected_clean)
# [1] 33848 17031

# Now making sure phenotype data from sce is in sce_uncorrected (clean)
colData(sce_uncorrected_clean) <- colData(sce)

# just checking for final Annotations
table(sce_uncorrected_clean$final_Annotations)
    # Astrocyte       Endo Excit.Thal Inhib.Thal      LHb.1      LHb.2      LHb.3 
    # 538         38       1800       7612        201        266        134 
    # LHb.4      LHb.5      LHb.6      LHb.7      MHb.1      MHb.2      MHb.3 
    # 477         83         39       1014        152        540         18 
    # Microglia      Oligo        OPC 
    # 145       2178       1796 

#### Prepping for plotting #####################################################
# creating plot directory
plot_dir <- here("plots", "99_paper_figs", "01_TSNEs", "Pre-Harmony")
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

##### PLOTTING TSNEs ###########################################################
# Pre-Harmonization colored by Samplel and faceted by Sample
pdf(here(plot_dir, "TSNE_uncorrected_by_Sample.pdf"), width = 14, height = 9)
plot1 <-  plotReducedDim(sce_uncorrected_clean, dimred = "TSNE") +
    geom_point(aes(color = sce_uncorrected_clean$Sample)) + 
    guides(color = guide_legend(title="Sample ID")) + 
  theme(legend.position = "none")

plot2 <-  plotReducedDim(sce_uncorrected_clean, dimred = "TSNE") +
  geom_point(aes(color = sce_uncorrected_clean$Sample)) + 
  guides(color = guide_legend(title="Sample ID")) +
  facet_wrap(~ sce_uncorrected_clean$Sample)

plot_grid(plot1, plot2)
dev.off()

# Pre-Harmonization colored by Samplel and faceted by Run
pdf(here(plot_dir, "TSNE_uncorrected_by_Run.pdf"), 
    width = 14, height = 9)
plot1 <- plotReducedDim(sce_uncorrected_clean, dimred = "TSNE") +
  geom_point(aes(color = sce_uncorrected_clean$Sample)) + 
  theme(legend.position = "none")

plot2 <- plotReducedDim(sce_uncorrected_clean, dimred = "TSNE") +
  geom_point(aes(color = sce_uncorrected_clean$Sample)) + 
  facet_wrap(~ sce_uncorrected_clean$Run) + 
  guides(color = guide_legend(title="Sample ID"))

plot_grid(plot1, plot2, rel_widths = c(1,2))
dev.off()

# for NeuN sorting and Non-NeuN sorting
plot_sorted <- plotReducedDim(sce_unc_sorted, dimred = "TSNE") +
  geom_point(aes(color = sce_unc_sorted$final_Annotations), alpha = 0.4) + 
  scale_colour_manual(values = cluster_colors) +
  facet_grid(sce_unc_sorted$NeuN ~ sce_unc_sorted$Sample) + 
  guides(color = guide_legend(title="Cell Type"))

plot_unsorted <-   plotReducedDim(sce_unc_unsorted, dimred = "TSNE") +
  geom_point(aes(color = sce_unc_unsorted$final_Annotations), alpha = 0.4) + 
  scale_colour_manual(values = cluster_colors) +
  facet_grid(sce_unc_unsorted$NeuN ~ sce_unc_unsorted$Sample) + 
  guides(color = guide_legend(title="Cell Type"))


pdf(here(plot_dir, "TSNE_harmony_by_finalAnno_splitbySampleAndSorting_PRE-HARMONY.pdf"), width = 13, height = 9)
plot_grid(
  plot_sorted,
  plot_unsorted,
  ncol = 1
)
dev.off()

pdf(here(plot_dir, "TSNE_harmony_by_final_Annotations_splitbyRun_PRE-HARMONY.pdf"), 
    width = 20, height = 8)
plot1 <- plotReducedDim(sce_uncorrected_clean, dimred = "TSNE") +
  geom_point(aes(color = sce_uncorrected_clean$final_Annotations)) + 
  scale_colour_manual(values = cluster_colors) +
  facet_wrap(~ sce_uncorrected_clean$Run) + 
  theme(legend.position = "none")

plot2 <- plotReducedDim(sce_uncorrected_clean, dimred = "TSNE") +
  geom_point(aes(color = sce_uncorrected_clean$final_Annotations)) + 
  scale_colour_manual(values = cluster_colors) +
  facet_wrap(~ sce_uncorrected_clean$final_Annotations) + 
  guides(color = guide_legend(title="Cell Type"))

plot_grid(plot1, plot2)
dev.off()

# Saving uncorrected sce object 
save(sce_uncorrected_clean, file = here("processed-data", "99_paper_figs", "sce_objects", 
           "official_final_uncorrected_sce.RDATA"))

# Done.

