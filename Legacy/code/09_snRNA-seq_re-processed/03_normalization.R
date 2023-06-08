library(SingleCellExperiment)
library(DropletUtils)
library(scater)
library(scran)
library(scry)
library(batchelor)
library(bluster)
library(BiocParallel)
library(jaffelab)
library(gridExtra)
library(here)
library(sessioninfo)
library(ggplot2)

load(here("processed-data","08_snRNA-seq_Erik", "01_qc.rda"), verbose = TRUE)

table(sce.all.hb$discard, sce.all.hb$discard_auto)
table(!sce.all.hb$discard &! sce.all.hb$discard_auto)
sce.all.hb <- sce.all.hb[,!sce.all.hb$discard &! sce.all.hb$discard_auto]
# Use `multiBatchNorm()` to compute log-normalized counts, matching the scaling across samples
sce.all.hb <- multiBatchNorm(sce.all.hb, batch=sce.all.hb$sample_short)

# Use the simple `modelGeneVar` - this makes more sense over `combineVar`, since the
#   cell composition is already known to be quite different (with NeuN selection)

# [1] 13076


### Dimensionality reduction ================================================================

# Run `fastMNN` (internally uses `multiBatchPCA`), taking top 100 (instead of default 50 PCs)
set.seed(109)

sce.all.hb <- devianceFeatureSelection(sce.all.hb,
                                   assay="counts", fam="binomial", sorted=F,
                                   batch=as.factor(sce.all.hb$sample_short))
    # This temp file just used for getting batch-corrected components (drops a variety of entries)

pdf(here("plots","09_snRNA-seq_re-processed", "binomial_deviance.pdf"))
plot(sort(rowData(sce.all.hb)$binomial_deviance, decreasing=T),
     type="l", xlab="ranked genes",
     ylab="binomial deviance", main="Feature Selection with Deviance: LC (n=3)")
abline(v=2000, lty=2, col="red")
dev.off()

Sys.time()
#"2022-06-16 12:15:58 EDT"
sce.all.hb <- nullResiduals(sce.all.hb, assay="counts", fam="binomial",  # default params
                        type="deviance")


hdgs.hb <- rownames(sce.all.hb)[order(rowData(sce.all.hb)$binomial_deviance, decreasing=T)][1:2000]

sce.all.hb <-  runPCA(sce.all.hb, exprs_values="binomial_deviance_residuals",
                  subset_row=hdgs.hb, ncomponents=100,
                  name="GLMPCA_approx",
                  BSPARAM=BiocSingular::IrlbaParam())

glmpca.mnn <- reducedMNN(reducedDim(sce.all.hb, "GLMPCA_approx"),
                         batch=as.factor(sce.all.hb$sample_short),
                         merge.order=c("Br1092", "Br1204", "Br1469", "Br1735", "Br5555", "Br5558", "Br5639")
                         )

sce.all.hb<- multiBatchNorm(sce.all.hb, batch=sce.all.hb$sample_short)
Sys.time()
    #[1] "2022-06-16 12:45:34 EDT"

# Store this
reducedDim(sce.all.hb, "GLMPCA_MNN") <- glmpca.mnn$corrected








pdf(file = here("plots","08_snRNA-seq_Erik", "GLMPCA_MNN_sample_id.pdf"), width = 9)
ggplot(
    data.frame(reducedDim(sce.all.hb, "GLMPCA_MNN")),
    aes(x = PC1, y = PC2, color = factor(sce.all.hb$sample_short))
) +
    geom_point() +
    labs(color = "Sample") +
    theme_bw()
dev.off()



pdf(file = here("plots","08_snRNA-seq_Erik", "uncorrected_sample_id.pdf"), width = 9)
ggplot(
    data.frame(reducedDim(sce.all.hb)),
    aes(x = PC1, y = PC2, color = factor(sce.all.hb$sample_short))
) +
    geom_point() +
    labs(color = "Sample") +
    theme_bw()
dev.off()


# Save into a new region-specific SCE object/flie
save(sce.all.hb,
     file=here("processed-data","09_snRNA-seq_re-processed","02_normalization.Rda"))


sessionInfo()
# R version 4.1.2 Patched (2021-11-04 r81138)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)
#
# Matrix products: default
# BLAS:   /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/lib/libRblas.so
# LAPACK: /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/lib/libRlapack.so
#
# locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
#  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
#  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
#  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
#
# attached base packages:
# [1] stats4    stats     graphics  grDevices datasets  utils     methods
# [8] base
#
# other attached packages:
#  [1] sessioninfo_1.2.2           here_1.0.1
#  [3] gridExtra_2.3               jaffelab_0.99.32
#  [5] rafalib_1.0.0               BiocParallel_1.28.3
#  [7] bluster_1.4.0               batchelor_1.10.0
#  [9] scry_1.6.0                  scran_1.22.1
# [11] scater_1.22.0               ggplot2_3.3.6
# [13] scuttle_1.4.0               DropletUtils_1.14.2
# [15] SingleCellExperiment_1.16.0 SummarizedExperiment_1.24.0
# [17] Biobase_2.54.0              GenomicRanges_1.46.1
# [19] GenomeInfoDb_1.30.1         IRanges_2.28.0
# [21] S4Vectors_0.32.4            BiocGenerics_0.40.0
# [23] MatrixGenerics_1.6.0        matrixStats_0.62.0
#
# loaded via a namespace (and not attached):
#  [1] segmented_1.4-1           fs_1.5.2
#  [3] bitops_1.0-7              RColorBrewer_1.1-3
#  [5] rprojroot_2.0.3           tools_4.1.2
#  [7] utf8_1.2.2                R6_2.5.1
#  [9] irlba_2.3.5               ResidualMatrix_1.4.0
# [11] HDF5Array_1.22.1          vipor_0.4.5
# [13] DBI_1.1.2                 colorspace_2.0-3
# [15] rhdf5filters_1.6.0        withr_2.5.0
# [17] tidyselect_1.1.2          compiler_4.1.2
# [19] cli_3.3.0                 BiocNeighbors_1.12.0
# [21] DelayedArray_0.20.0       labeling_0.4.2
# [23] scales_1.2.0              digest_0.6.29
# [25] R.utils_2.11.0            XVector_0.34.0
# [27] pkgconfig_2.0.3           sparseMatrixStats_1.6.0
# [29] limma_3.50.3              rlang_1.0.2
# [31] DelayedMatrixStats_1.16.0 farver_2.1.0
# [33] generics_0.1.2            dplyr_1.0.9
# [35] R.oo_1.25.0               RCurl_1.98-1.7
# [37] magrittr_2.0.3            BiocSingular_1.10.0
# [39] GenomeInfoDbData_1.2.7    Matrix_1.4-1
# [41] Rcpp_1.0.8.3              ggbeeswarm_0.6.0
# [43] munsell_0.5.0             Rhdf5lib_1.16.0
# [45] fansi_1.0.3               viridis_0.6.2
# [47] lifecycle_1.0.1           R.methodsS3_1.8.2
# [49] edgeR_3.36.0              MASS_7.3-56
# [51] zlibbioc_1.40.0           rhdf5_2.38.1
# [53] grid_4.1.2                parallel_4.1.2
# [55] ggrepel_0.9.1             dqrng_0.3.0
# [57] crayon_1.5.1              lattice_0.20-45
# [59] beachmat_2.10.0           splines_4.1.2
# [61] locfit_1.5-9.5            metapod_1.2.0
# [63] pillar_1.7.0              igraph_1.3.2
# [65] ScaledMatrix_1.2.0        glue_1.6.2
# [67] vctrs_0.4.1               gtable_0.3.0
# [69] purrr_0.3.4               assertthat_0.2.1
# [71] rsvd_1.0.5                googledrive_2.0.0
# [73] viridisLite_0.4.0         gargle_1.2.0
# [75] tibble_3.1.7              beeswarm_0.4.0
# [77] cluster_2.1.3             statmod_1.4.36
# [79] ellipsis_0.3.2
