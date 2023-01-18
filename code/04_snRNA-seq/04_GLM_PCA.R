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

# loading post qc sce object 
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
           "sce_hb_postQC.Rdata"))
sce <- sce_hb_postQC
rm(sce_hb_postQC)

# Deviance featuring selection
sce <- devianceFeatureSelection(sce,
        assay = "counts", fam = "binomial", sorted = F,
        batch = as.factor(sce$Sample))

# Checking outputs for deviance ft selection function
pdf(here("plots", "04_snRNA-seq", "04_GLM_PCA_plots", "binomial_deviance.pdf"))
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

# Running PCA 
sce_uncorrected <- runPCA(sce,
                          exprs_values = "binomial_deviance_residuals",
                          subset_row = hdgs.hb, ncomponents = 100,
                          name = "GLMPCA_approx",
                          BSPARAM = BiocSingular::IrlbaParam()
)

# Plotting using GLM-PCA functions 
# We should color by batches but we must track that down
pdf(here("plots", "04_snRNA-seq", "04_GLM_PCA_plots", "main_GLM_pca_plot.pdf"))
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "Sample") # plot PC1&2, facet_wrap by sample
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "Sample", ncomponents = c(2,3)) # looks batchy
  plotReducedDim(sce_uncorrected, dimred = "GLMPCA_approx", colour_by = "Sample", ncomponents = 5)
dev.off()

# Running TSNE and UMAP
sce_uncorrected <- runTSNE(sce_uncorrected, dimred = "GLMPCA_approx")
sce_uncorrected <- runUMAP(sce_uncorrected, dimred = "GLMPCA_approx")

# Plotting TSNE
pdf(here("plots", "04_snRNA-seq", "04_GLM_PCA_plots", "TSNE_uncorrected_plot.pdf"))
  plotReducedDim(sce_uncorrected, dimred = "TSNE", colour_by = "Sample") 
dev.off()

# Plotting UMAP
# Plotting TSNE
pdf(here("plots", "04_snRNA-seq", "04_GLM_PCA_plots", "UMAP_uncorrected_plot.pdf"))
  plotReducedDim(sce_uncorrected, dimred = "UMAP", colour_by = "Sample") 
dev.off()






## Save uncorrected sce object post pca 
save(sce_uncorrected, file = here("processed-data", "04_snRNA-seq", "sce_objects", 
                                  "sce_uncorrected_PCA.Rdata"))

# Saving as HDF5 for later clustering
saveHDF5SummarizedExperiment(sce_uncorrected, dir = here("processed-data", "04_snRNA-seq", 
                            "sce_objects", "sce_uncorrected"))




## Q for Louise: Should we color by cell type? And why does our sce object not have 
# cell-type? Is that something we compute on our own or did we lose this in our various
# saves?

# Steps: 1) Normalization, 2) Feature Selection (get rid of features with low variance),
# 3) Dimensionality Reduction, 4) Clustering, 5) Marker gene detection,
# 6) Cell Type Annotation


# TO DO:
# ** Push with Collabo
# Slack Write up 
# 1) Track down batch info in code to color by.
# 2) Plot everything again by batch data but also qc metrics.
# 3) Share on Slack 
# 4) Run harmony by sample.
# 5) 1/18/23 (Louise): Work on doublet detection for sn qc.

