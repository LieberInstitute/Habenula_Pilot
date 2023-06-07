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

## Use previous norm factors to scale the raw library sizes
rse_gene_filt_scaled <- calcNormFactors(rse_gene_filt)
rse_gene_filt_scaled$samples$library_size <- factors$library_size
rse_gene_filt_scaled$samples$norm.factors <- factors$norm.factors

## Transform counts to log2(CPM)
## Estimate mean-variance relationship for each gene
vGene <- voom(rse_gene_filt_scaled, design = model, plot = TRUE)

## Fit linear model for each gene
fitGene <- lmFit(vGene)

## Empirical Bayesian calculation to obtain the significant genes:
## compute moderated F and t-statistics, and log-odds of DE
eBGene <- eBayes(fitGene)

## Plot average log expression vs logFC
limma::plotMA(eBGene,
    coef = "PrimaryDxSchizo",
    xlab = "Mean of normalized counts",
    ylab = "logFC"
)

## Plot -log(p-value) vs logFC
volcanoplot(eBGene, coef = "PrimaryDxSchizo")

## Top-ranked genes for Substance (cases vs ctrls)
top_genes <- topTable(eBGene, coef = "PrimaryDxSchizo", p.value = 1, number = nrow(rse_gene_filt), sort.by = "none")

## Histogram of adjusted p values
hist(top_genes$adj.P.Val, xlab = "FDR", main = "")

###############################################################################



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
