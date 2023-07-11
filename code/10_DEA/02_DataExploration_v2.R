library("here")
library("SummarizedExperiment")
library("PCAtools")
library("stringr")
library("ggplot2")
library("cowplot")
library("GGally")
library("ComplexHeatmap")
library("circlize")
library("sessioninfo")

output_path <- here("plots", "10_DEA", "02_DataExploration")



############################# Load rse gene object ############################

load(
    here(
        "processed-data",
        "rse_objects",
        "rse_gene_Habenula_Pilot.rda"
    ),
)

lobstr::obj_size(rse_gene)
# 29.92 MB

class(rse_gene)
# [1] "RangedSummarizedExperiment"
# attr(,"package")
# [1] "SummarizedExperiment"

dim(rse_gene)
# [1] 22756    68

###############################################################################




################## PCA and heatmap with colData() variables ###################

## NOTE: This section is done with filtered and normalized counts

## PCA
pca_df <- pca(assays(rse_gene)$logcounts, metadata = colData(rse_gene))

## Plot correlation between PCs and variables (with and without qSVs)
pdf(paste(output_path, "/", "Corr_PCA-Vars-qSVs.pdf", sep = ""), width = 10, height = 10)
eigencorplot(pca_df, metavars = c("PrimaryDx", "AgeDeath", "Flowcell", "mitoRate", "rRNA_rate", "totalAssignedGene", "RIN", "abs_ERCCsumLogErr", "snpPC1", "snpPC2", "snpPC3", "snpPC4", "snpPC5", "tot.Hb", "tot.Thal", "qSV1", "qSV2", "qSV3", "qSV4", "qSV5"))
dev.off()

pdf(paste(output_path, "/", "Corr_PCA-Vars-noqSVs.pdf", sep = ""), width = 10, height = 10)
eigencorplot(pca_df, metavars = c("PrimaryDx", "AgeDeath", "Flowcell", "mitoRate", "rRNA_rate", "totalAssignedGene", "RIN", "abs_ERCCsumLogErr", "snpPC1", "snpPC2", "snpPC3", "snpPC4", "snpPC5", "tot.Hb", "tot.Thal"))
dev.off()

###############################################################################



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
