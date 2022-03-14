library(here)
library(sessioninfo)
library(batchelor)
library(scran)

load(here("processed-data","08_snRNA-seq_Erik", "20220301_human_hb_processing.rda"), verbose = TRUE)


# Use `multiBatchNorm()` to compute log-normalized counts, matching the scaling across samples
sce.all.hb <- multiBatchNorm(sce.all.hb, batch=sce.all.hb$sample_short)

# Use the simple `modelGeneVar` - this makes more sense over `combineVar`, since the
#   cell composition is already known to be quite different (with NeuN selection)
geneVar <- modelGeneVar(sce.all.hb)
chosen.hvgs.hb <- geneVar$bio > 0
sum(chosen.hvgs.hb)
# [1] 13076


### Dimensionality reduction ================================================================

# Run `fastMNN` (internally uses `multiBatchPCA`), taking top 100 (instead of default 50 PCs)
set.seed(109)
mnn.hold <-  fastMNN(sce.all.hb, batch=sce.all.hb$sample_short,
                     merge.order=c("Br1092","Br1204","Br1469","Br1735","Br5555","Br5558","Br5639"),
                     subset.row=chosen.hvgs.hb, d=100,
                     correct.all=TRUE, get.variance=TRUE,
                     BSPARAM=BiocSingular::IrlbaParam())
    # This temp file just used for getting batch-corrected components (drops a variety of entries)

table(colnames(mnn.hold) == colnames(sce.all.hb))  # all TRUE
table(mnn.hold$batch == sce.all.hb$sample_short) # all TRUE

# Add them to the SCE, as well as the metadata (though the latter might not be so usefl)
reducedDim(sce.all.hb, "PCA_corrected") <- reducedDim(mnn.hold, "corrected") # 100 components
metadata(sce.all.hb) <- metadata(mnn.hold)

# Save into a new region-specific SCE object/flie
save(sce.all.hb, chosen.hvgs.hb,
     file="/dcl02/lieber/ajaffe/Roche_Habenula/processed-data/09_snRNA-seq_re-processed/02_normalization.Rda")
