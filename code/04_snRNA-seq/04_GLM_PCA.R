# January 17, 2023 - Bukola Ajanaku
# Creating pca plots for filtered (removed empty droplets), quality controlled 
# (dropped high mito, low library size, low detected features, and any genes 
# with 0 counts across samples) sce object.
# Based on:
# 1) https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/code/03_build_sce/04_GLM_PCA.R
# 2) https://www.stephaniehicks.com/biocdemo/articles/Demo.html
# qrsh -l mem_free=50G,h_vmem=50G

library("SingleCellExperiment")
library("scran")
library("scater")
library("scry")
library("here")
library("sessioninfo")
library("mbkmeans")
library("HDF5Array")
library("ggplot2")
library("dplyr")

# loading post qc sce object 
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
           "sce_hb_postQC.Rdata"))
sce <- sce_hb_postQC
rm(sce_hb_postQC)

####### Adding new metric: Run Data (before subsetting) ########################
# Adding Run data 
colData(sce)$Run <- NA
colData(sce)[colData(sce)$Sample == "Br1469", "Run"] <- 1
colData(sce)[colData(sce)$Sample %in% c("Br5558", "Br1204"), "Run"] <- 2
colData(sce)[colData(sce)$Sample %in% c("Br1092", "Br1735", "Br5555", "Br5639"), "Run"] <- 3
################################################################################

# Deviance featuring selection
sce <- devianceFeatureSelection(sce,
        assay = "counts", fam = "binomial", sorted = F,
        batch = as.factor(sce$Sample))

# Checking outputs for deviance ft selection function
pdf(here("plots", "04_snRNA-seq", "04_GLM_PCA_plots", "binomial_deviance_uncorrected.pdf"))
  plot(sort(rowData(sce)$binomial_deviance, decreasing = T),
       type = "l", xlab = "ranked genes",
       ylab = "binomial deviance", main = "Feature Selection with Deviance"
  )
  abline(v = 2000, lty = 2, col = "red")
dev.off()

# Taking the GLM_PCA approach
sce <- nullResiduals(sce,
                     assay = "counts", fam = "binomial", # default params
                     type = "deviance")

# Selects for the top 2000 most variable genes
hdgs.hb <- rownames(sce)[order(rowData(sce)$binomial_deviance, decreasing = T)][1:2000]

# RUNNING PCA --
sce_uncorrected <- runPCA(sce,
                          exprs_values = "binomial_deviance_residuals",
                          subset_row = hdgs.hb, ncomponents = 100,
                          name = "GLMPCA_approx",
                          BSPARAM = BiocSingular::IrlbaParam()
)

# Plotting using GLM-PCA functions [color by Sample]
pdf(here("plots", "04_snRNA-seq", "04_GLM_PCA_plots", "GLM_uncorrected_plot_by_Sample.pdf"))
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "Sample") 
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "Sample") + facet_wrap(~ sce_uncorrected$Sample)
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "Sample", ncomponents = c(2,3)) # looks batchy
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "Sample", ncomponents = c(2,3)) + 
    facet_wrap(~ sce_uncorrected$Sample)
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "Sample", ncomponents = 5)
dev.off()

# Plotting using GLM-PCA functions [color by Run]
pdf(here("plots", "04_snRNA-seq", "04_GLM_PCA_plots", "GLM_uncorrected_plot_by_Run.pdf"))
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "Run") 
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "Run") + facet_wrap(~ sce_uncorrected$Sample)
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "Run", ncomponents = c(2,3)) # looks batchy
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "Run", ncomponents = c(2,3)) + 
    facet_wrap(~ sce_uncorrected$Sample)
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "Run", ncomponents = 5)
dev.off()

# Plotting by continuous data
pdf(here("plots", "04_snRNA-seq", "04_GLM_PCA_plots", "GLM_uncorrected_continous_metrics.pdf"))
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "sum") 
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "sum")  + facet_wrap(~ sce_uncorrected$Sample)
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "sum")  + facet_wrap(~ sce_uncorrected$Run)
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "detected") 
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "detected")  + facet_wrap(~ sce_uncorrected$Sample)
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "detected")  + facet_wrap(~ sce_uncorrected$Run)
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "subsets_Mito_percent") 
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "subsets_Mito_percent")  + facet_wrap(~ sce_uncorrected$Sample)
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "subsets_Mito_percent")  + facet_wrap(~ sce_uncorrected$Run)
dev.off()

# RUNNING TSNE and UMAP --
sce_uncorrected <- runTSNE(sce_uncorrected, dimred = "GLMPCA_approx")
sce_uncorrected <- runUMAP(sce_uncorrected, dimred = "GLMPCA_approx")

# Plotting TSNE 
pdf(here("plots", "04_snRNA-seq", "04_GLM_PCA_plots", "TSNE_uncorrected_plot.pdf"))
  plotReducedDim(sce_uncorrected, dimred = "TSNE", colour_by = "Sample") 
  plotReducedDim(sce_uncorrected, dimred = "TSNE", colour_by = "Sample") + facet_wrap(~ sce_uncorrected$Run)
  plotReducedDim(sce_uncorrected, dimred = "TSNE", colour_by = "Run") 
  plotReducedDim(sce_uncorrected, dimred = "TSNE", colour_by = "Run") + facet_wrap(~ sce_uncorrected$Sample)
dev.off()

# Plotting UMAP
pdf(here("plots", "04_snRNA-seq", "04_GLM_PCA_plots", "UMAP_uncorrected_plot.pdf"))
  plotReducedDim(sce_uncorrected, dimred = "UMAP", colour_by = "Sample")
  plotReducedDim(sce_uncorrected, dimred = "UMAP", colour_by = "Sample") + facet_wrap(~ sce_uncorrected$Run)
  plotReducedDim(sce_uncorrected, dimred = "UMAP", colour_by = "Run") 
  plotReducedDim(sce_uncorrected, dimred = "UMAP", colour_by = "Run") + facet_wrap(~ sce_uncorrected$Sample)
dev.off()

# Plotting by continuous data, TSNE
pdf(here("plots", "04_snRNA-seq", "04_GLM_PCA_plots", "TSNE_uncorrected_plot_continous_mets.pdf"))
  plotReducedDim(sce_uncorrected, dimred = "TSNE", colour_by = "sum") 
  plotReducedDim(sce_uncorrected, dimred = "TSNE", colour_by = "sum") + facet_wrap(~ sce_uncorrected$Sample)
  plotReducedDim(sce_uncorrected, dimred = "TSNE", colour_by = "sum") + facet_wrap(~ sce_uncorrected$Run)
  plotReducedDim(sce_uncorrected, dimred = "TSNE", colour_by = "detected") 
  plotReducedDim(sce_uncorrected, dimred = "TSNE", colour_by = "detected") + facet_wrap(~ sce_uncorrected$Sample)
  plotReducedDim(sce_uncorrected, dimred = "TSNE", colour_by = "detected") + facet_wrap(~ sce_uncorrected$Run)
  plotReducedDim(sce_uncorrected, dimred = "TSNE", colour_by = "detected") 
  plotReducedDim(sce_uncorrected, dimred = "TSNE", colour_by = "subsets_Mito_percent") + facet_wrap(~ sce_uncorrected$Sample)
  plotReducedDim(sce_uncorrected, dimred = "TSNE", colour_by = "subsets_Mito_percent") + facet_wrap(~ sce_uncorrected$Run)
dev.off()

# Plotting by continous data, UMAP
pdf(here("plots", "04_snRNA-seq", "04_GLM_PCA_plots", "UMAP_uncorrected_plot_continous_mets.pdf"))
  plotReducedDim(sce_uncorrected, dimred = "UMAP", colour_by = "sum") 
  plotReducedDim(sce_uncorrected, dimred = "UMAP", colour_by = "sum") + facet_wrap(~ sce_uncorrected$Sample)
  plotReducedDim(sce_uncorrected, dimred = "UMAP", colour_by = "sum") + facet_wrap(~ sce_uncorrected$Run)
  plotReducedDim(sce_uncorrected, dimred = "UMAP", colour_by = "detected") 
  plotReducedDim(sce_uncorrected, dimred = "UMAP", colour_by = "detected") + facet_wrap(~ sce_uncorrected$Sample)
  plotReducedDim(sce_uncorrected, dimred = "UMAP", colour_by = "detected") + facet_wrap(~ sce_uncorrected$Run)
  plotReducedDim(sce_uncorrected, dimred = "UMAP", colour_by = "detected") 
  plotReducedDim(sce_uncorrected, dimred = "UMAP", colour_by = "subsets_Mito_percent") + facet_wrap(~ sce_uncorrected$Sample)
  plotReducedDim(sce_uncorrected, dimred = "UMAP", colour_by = "subsets_Mito_percent") + facet_wrap(~ sce_uncorrected$Run)
dev.off()


## Save uncorrected sce object post pca 
save(sce_uncorrected, file = here("processed-data", "04_snRNA-seq", "sce_objects", 
                                  "sce_uncorrected_PCA.Rdata"))

# Saving as HDF5 for later clustering
saveHDF5SummarizedExperiment(sce_uncorrected, dir = here("processed-data", "04_snRNA-seq", 
                            "sce_objects", "sce_uncorrected"))




# TO DO:
# 1) Track down batch info in code to color by.
# 2) Plot everything again by batch data but also qc metrics.
# 3) Share on Slack 
# 4) Run harmony by sample.
# 5) 1/18/23 (Louise): Work on doublet detection for sn qc.

