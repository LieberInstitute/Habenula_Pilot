# January 18, 2023 - Bukola Ajanaku
# Normalizing sce object by different metrics and then plotting PCs, TSNEs, and
# UMAPs.
# Based on:
# 1) https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/code/03_build_sce/05_harmony_correction.R
# 2) https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/code/03_build_sce/05.5_harmony_correction_plots.R
# qrsh -l mem_free=50G,h_vmem=50G

# Loading relevant libraries
library("SingleCellExperiment")
library("harmony")
library("scater")
library("here")
library("sessioninfo")
library("HDF5Array")
library("viridis")
library("ggplot2")


# Loading sce_uncorrected
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
          "sce_uncorrected_PCA.Rdata"))

####### RUNNING HARMONY ########################################################
# Note: We used GLM PCA not normal PCA but harmony function searches for something
# named "PCA". Thus, we must submit GLM under the name PCA. After running harmony,
# we can do away with our fake PCA:

# Submitting GLM PCA under the name PCA:
reducedDim(sce_uncorrected, "PCA") <- reducedDim(sce_uncorrected, "GLMPCA_approx")

# Running harmony
sce_corrbySamp <- RunHarmony(sce_uncorrected, group.by.vars = "Sample", verbose = TRUE)
# sce_corrbyRun <- RunHarmony(sce_uncorrected, group.by.vars = "Run", verbose = TRUE)
    # Error in harmonyObj$init_cluster_cpp(0) : 
    # element-wise multiplication: incompatible matrix dimensions: 100x3 and 100x1

# Removing our redudant fake PCA name and keeping it under "GLMPCA_approx"
reducedDim(sce_uncorrected, "PCA") <- NULL

####### TSNE & UMAP ############################################################
set.seed(777)

## sce corrected by Sample
sce_corrbySamp <- runTSNE(sce_corrbySamp, dimred = "HARMONY")
sce_corrbySamp <- runUMAP(sce_corrbySamp, dimred = "HARMONY")

## sce corrected by Run (not possible because uneven)
# sce_corrbyRun <- runTSNE(sce_corrbyRun, dimred = "HARMONY")
# sce_corrbyRun <- runUMAP(sce_corrbyRun, dimred = "HARMONY")

####### PLOTTING ###############################################################
# GLM by Sample
pdf(here("plots", "04_snRNA-seq", "05_GLM_Harmony_plots", "GLM_harmony_plot_by_Sample.pdf"))
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "Sample") 
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "Sample") + 
    facet_wrap(~ sce_corrbySamp$Sample)
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "Sample", ncomponents = c(2,3)) # looks batchy
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "Sample", ncomponents = c(2,3)) + 
    facet_wrap(~ sce_corrbySamp$Sample)
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "Sample", ncomponents = 5)
dev.off()

# Plotting using GLM-PCA functions [color by Run]
pdf(here("plots", "04_snRNA-seq", "05_GLM_Harmony_plots", "GLM_harmony_plot_by_Run.pdf"))
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "Run") 
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "Run") + 
    facet_wrap(~ sce_corrbySamp$Sample)              
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "Run", ncomponents = c(2,3)) 
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "Run", ncomponents = c(2,3)) + 
    facet_wrap(~ sce_corrbySamp$Sample)
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "Run", ncomponents = 5)
dev.off()

# Plotting using GLM-PCA functions [color by Erik Cluster]
pdf(here("plots", "04_snRNA-seq", "05_GLM_Harmony_plots", "GLM_harmony_plot_by_ct_Erik.pdf"))
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "ct_Erik") + labs(caption = "887 NAs")
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "ct_Erik") + 
    facet_wrap(~ sce_corrbySamp$Sample)  + labs(caption = "887 NAs")            
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "ct_Erik", ncomponents = c(2,3)) + 
    labs(caption = "887 NAs")
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "ct_Erik", ncomponents = c(2,3)) + 
    facet_wrap(~ sce_corrbySamp$Sample) + labs(caption = "887 NAs")
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "ct_Erik", ncomponents = 5) + 
    labs(caption = "887 NAs")
dev.off()

# GLM by Continous Metrics 
pdf(here("plots", "04_snRNA-seq", "05_GLM_Harmony_plots", "GLM_harmony_continous_metrics.pdf"))
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "sum") +
    scale_color_viridis(option = "G", begin = 0.19) 
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "sum")  + 
    facet_wrap(~ sce_corrbySamp$Sample) +
    scale_color_viridis(option = "G", begin = 0.19) 
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "detected") +
    scale_color_viridis(option = "G", begin = 0.19) 
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "detected")  + 
    facet_wrap(~ sce_corrbySamp$Sample) +
    scale_color_viridis(option = "G", begin = 0.19) 
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "subsets_Mito_percent") +
    scale_color_viridis(option = "G", begin = 0.19) 
  plotReducedDim(sce_corrbySamp, dimred = "HARMONY", colour_by = "subsets_Mito_percent") + 
    facet_wrap(~ sce_corrbySamp$Sample) +
    scale_color_viridis(option = "G", begin = 0.19) 
dev.off()

# TSNE by Sample
pdf(here("plots", "04_snRNA-seq", "05_GLM_Harmony_plots", "TSNE_harmony_plot.pdf"))
  plotReducedDim(sce_corrbySamp, dimred = "TSNE", colour_by = "Sample") 
  plotReducedDim(sce_corrbySamp, dimred = "TSNE", colour_by = "Sample") + facet_wrap(~ sce_corrbySamp$Sample)
  plotReducedDim(sce_corrbySamp, dimred = "TSNE", colour_by = "Run") 
  plotReducedDim(sce_corrbySamp, dimred = "TSNE", colour_by = "Run") + facet_wrap(~ sce_corrbySamp$Sample)
dev.off()

# TSNE by continuous data
pdf(here("plots", "04_snRNA-seq", "05_GLM_Harmony_plots", "TSNE_harmony_plot_continous_mets.pdf"))
  plotReducedDim(sce_corrbySamp, dimred = "TSNE", colour_by = "sum") 
  plotReducedDim(sce_corrbySamp, dimred = "TSNE", colour_by = "sum") + facet_wrap(~ sce_corrbySamp$Sample)
  plotReducedDim(sce_corrbySamp, dimred = "TSNE", colour_by = "detected") 
  plotReducedDim(sce_corrbySamp, dimred = "TSNE", colour_by = "detected") + facet_wrap(~ sce_corrbySamp$Sample)
  plotReducedDim(sce_corrbySamp, dimred = "TSNE", colour_by = "subsets_Mito_percent") 
  plotReducedDim(sce_corrbySamp, dimred = "TSNE", colour_by = "subsets_Mito_percent") + facet_wrap(~ sce_corrbySamp$Sample)
dev.off()

# UMAP by Sample 
pdf(here("plots", "04_snRNA-seq", "05_GLM_Harmony_plots", "UMAP_harmony_plot.pdf"))
  plotReducedDim(sce_corrbySamp, dimred = "UMAP", colour_by = "Sample")
  plotReducedDim(sce_corrbySamp, dimred = "UMAP", colour_by = "Sample") + facet_wrap(~ sce_corrbySamp$Sample)
  plotReducedDim(sce_corrbySamp, dimred = "UMAP", colour_by = "Run") 
  plotReducedDim(sce_corrbySamp, dimred = "UMAP", colour_by = "Run") + facet_wrap(~ sce_corrbySamp$Sample)
dev.off()

# UMAP by continuous data 
pdf(here("plots", "04_snRNA-seq", "05_GLM_Harmony_plots", "UMAP_harmony_plot_continous_mets.pdf"))
  plotReducedDim(sce_corrbySamp, dimred = "UMAP", colour_by = "sum") 
  plotReducedDim(sce_corrbySamp, dimred = "UMAP", colour_by = "sum") + facet_wrap(~ sce_corrbySamp$Sample)
  plotReducedDim(sce_corrbySamp, dimred = "UMAP", colour_by = "detected") 
  plotReducedDim(sce_corrbySamp, dimred = "UMAP", colour_by = "detected") + facet_wrap(~ sce_corrbySamp$Sample)
  plotReducedDim(sce_corrbySamp, dimred = "UMAP", colour_by = "subsets_Mito_percent") 
  plotReducedDim(sce_corrbySamp, dimred = "UMAP", colour_by = "subsets_Mito_percent") + facet_wrap(~ sce_corrbySamp$Sample)
dev.off()


####### saving ###############################################################
## Saving harmonized (by Sample) sce object 
save(sce_corrbySamp, file = here("processed-data", "04_snRNA-seq", "sce_objects", 
                                  "sce_harmony_by_Samp.Rdata"))

# Saving harmonized (by Sample) sce object (Saved as HDF5 for later clustering)
saveHDF5SummarizedExperiment(sce_corrbySamp, dir = here("processed-data", "04_snRNA-seq", 
                                  "sce_objects", "sce_harmony_by_Samp"))




