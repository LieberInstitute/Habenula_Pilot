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
load(here("processed-data","08_snRNA-seq_Erik", "20220301_human_hb_processing.rda"), verbose = TRUE)


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
reducedDim(sce.lc, "GLMPCA_MNN") <- glmpca.mnn$corrected


table(colnames(mnn.hold) == colnames(sce.all.hb))  # all TRUE
table(mnn.hold$batch == sce.all.hb$sample_short) # all TRUE

# Add them to the SCE, as well as the metadata (though the latter might not be so usefl)
reducedDim(sce.all.hb, "PCA_corrected") <- reducedDim(mnn.hold, "corrected") # 100 components
metadata(sce.all.hb) <- metadata(mnn.hold)

# Save into a new region-specific SCE object/flie
save(sce.all.hb,
     file="/dcl02/lieber/ajaffe/Roche_Habenula/processed-data/09_snRNA-seq_re-processed/02_normalization.Rda")
