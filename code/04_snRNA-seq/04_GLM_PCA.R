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

# loading post qc sce object 
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
           "sce_hb_postQC.Rdata"))
sce <- sce_hb_postQC
rm(sce_hb_postQC)

# Normalization?
mbk <- mbkmeans(sce, whichAssay = "counts", reduceMethod = NA,
                clusters=10, batch_size = 500)
sce$mbk10 <- paste0("mbk", mbk$Clusters)
table(mbk$Clusters)

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

# Running PCA 
sce <- scater::runPCA(sce, ncomponents = 50,
                      exprs_values = "binomial_deviance_residuals",
                      scale = TRUE, name = "GLMPCA_approx",
                      BSPARAM = BiocSingular::RandomParam())

# Plotting using GLM-PCA functions
pdf(here("plots", "04_snRNA-seq", "04_GLM_PCA_plots", "GLM_pca_plot.pdf"))
  plotReducedDim(sce, dimred = "GLMPCA_approx", colour_by = "Sample")
  plotReducedDim(sce, dimred = "GLMPCA_approx", colour_by = "mbk10")
dev.off()

# Normal PCA
# had to re-load sce 
sce <- computeSumFactors(sce, cluster = mbk$Clusters, min.mean = 0.1)
sce <- logNormCounts(sce)

sce <- scater::runPCA(sce, ncomponents = 50,
                      ntop = 1000,
                      scale = TRUE,
                      BSPARAM = BiocSingular::RandomParam())

sce$mbk10 <- paste0("mbk", mbk$Clusters)

pdf(here("plots", "04_snRNA-seq", "04_GLM_PCA_plots", "normal_pca_plot.pdf"))
  plotPCA(sce, colour_by = "Sample")
  plotPCA(sce, colour_by = "mbk10")
dev.off()

## Q for Louise: Should we color by cell type? And why does our sce object not have 
# cell-type? Is that something we compute on our own or did we lose this in our various
# saves?

# Steps: 1) Normalization, 2) Feature Selection (get rid of features with low variance),
# 3) Dimensionality Reduction, 4) Clustering, 5) Marker gene detection,
# 6) Cell Type Annotation








## 