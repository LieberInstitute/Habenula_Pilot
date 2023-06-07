library(here)
library(SummarizedExperiment)
library(edgeR)
library(limma)
library(sessioninfo)



###################### Load rse and rse filtered objects ######################

load(
    here(
        "preprocessed_data",
        "count_data_bukola",
        "rse_gene_Roche_Habenula_qcAndAnnotated_n69.Rdata"
    ),
    verbose = TRUE
)
# Loading objects:
#   rse_gene

colData(rse_gene)$library_size <- apply(assay(rse_gene), 2, sum)
colData(rse_gene)$log10_library_size <- log10(colData(rse_gene)$library_size)
colData(rse_gene)$detected_num_genes <- apply(assay(rse_gene), 2, function(x) {
    length(x[which(x > 0)])
})
colData(rse_gene)$abs_ERCCsumLogErr <- abs(colData(rse_gene)$ERCCsumLogErr)

load(
    here(
        "processed-data",
        "10_DEA",
        "rse_gene_filt.Rdata"
    ),
    verbose = TRUE
)
# Loading objects:
#   rse_gene_filt

###############################################################################



####################### Differential Expression Analysis ######################

## Extract calcNormFactors for all samples
norm_factors <- calcNormFactors(rse_gene, method = "TMM")

samples_factors <- data.frame(
    BrNum = norm_factors$samples$BrNum,
    norm.factors = norm_factors$samples$norm.factors,
    library_size = norm_factors$samples$library_size
)

## Previous lib sizes of each sample
match_samples <- match(rse_gene_filt$BrNum, samples_factors$BrNum)
factors <- samples_factors[match_samples, ]

## Formula and model
formula <- ~ PrimaryDx + AgeDeath + Flowcell + mitoRate + rRNA_rate + totalAssignedGene + RIN + abs_ERCCsumLogErr + library_size
model <- model.matrix(formula, data = colData(rse_gene_filt))
######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
