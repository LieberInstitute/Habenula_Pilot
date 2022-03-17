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

pdf(here("processed-data","09_snRNA-seq_re-processed", "binomial_deviance.pdf"))
plot(sort(rowData(sce.all.hb)$binomial_deviance, decreasing=T),
     type="l", xlab="ranked genes",
     ylab="binomial deviance", main="Feature Selection with Deviance: LC (n=3)")
abline(v=2000, lty=2, col="red")
dev.off()

Sys.time()

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
    #[1] "2021-12-31 17:00:37 EST"

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
