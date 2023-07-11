library("here")
library("SummarizedExperiment")
library("edgeR")
library("limma")
library("dplyr")
library("EnhancedVolcano")
library("sessioninfo")

output_path <- here("plots", "10_DEA", "04_DEA")



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
rse_gene_raw <- rse_gene

load(
    here(
        "processed-data",
        "rse_objects",
        "rse_gene_Habenula_Pilot.rda"
    ),
    verbose = TRUE
)
# Loading objects:
#   rse_gene

###############################################################################



####################E### Functions for DEA and plotting #######################

norm_factors <- calcNormFactors(rse_gene_raw, method = "TMM")
samples_factors <- data.frame(
    RNum = norm_factors$samples$RNum,
    norm.factors = norm_factors$samples$norm.factors,
    lib.size = norm_factors$samples$lib.size
)


## Function to do all the DEA analysis with limma
DE_analysis <- function(rse_gene, formula, coef, model_name) {
    match_samples <- match(rse_gene$RNum, samples_factors$RNum)
    stopifnot(all(!is.na(match_samples)))
    factors <- samples_factors[match_samples, ]

    pdf(file = paste0(output_path, "/DEAplots_", model_name, ".pdf"))
    par(mfrow = c(2, 2))

    ## Model matrix
    model <- model.matrix(formula, data = colData(rse_gene))

    ## Use previous norm factors to scale the raw library sizes
    RSE_scaled <- calcNormFactors(rse_gene)
    RSE_scaled$samples$lib.size <- factors$lib.size
    RSE_scaled$samples$norm.factors <- factors$norm.factors

    ## Transform counts to log2(CPM): estimate mean-variance relationship for
    ## each gene
    vGene <- voom(RSE_scaled, design = model, plot = TRUE)

    ## Fit linear model for each gene
    fitGene <- lmFit(vGene)

    ## Empirical Bayesian calculation to obtain our significant genes: compute
    ## moderated F and t-statistics, and log-odds of DE
    eBGene <- eBayes(fitGene)

    ## Plot average log expression vs logFC
    limma::plotMA(eBGene,
        coef = coef, xlab = "Mean of normalized counts",
        ylab = "logFC"
    )

    ## Plot -log(p-value) vs logFC
    volcanoplot(eBGene, coef = coef)

    ## Select top-ranked genes
    top_genes <- topTable(eBGene, coef = coef, p.value = 1, number = nrow(rse_gene), sort.by = "none")

    ## Histogram of adjusted p values
    hist(top_genes$adj.P.Val, xlab = "FDR", main = "")

    dev.off()

    return(top_genes)
}

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

## Top-ranked genes for PrimaryDx (Schizo vs Control)
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
