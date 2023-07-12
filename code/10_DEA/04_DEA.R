library("here")
library("SummarizedExperiment")
library("edgeR")
library("limma")
library("dplyr")
library("EnhancedVolcano")
library("sessioninfo")

out_plot <- here("plots", "10_DEA", "04_DEA")
out_data <- here("processed-data", "10_DEA", "04_DEA")



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



#################### E### Functions for DEA and plotting #######################

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

    pdf(file = paste0(out_plot, "/DEAplots_", model_name, ".pdf"))
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

    write.table(top_genes, file = paste0(out_data, "/DEA_AllGenes_", model_name, ".tsv"), sep = "\t", quote = FALSE)

    return(top_genes)
}

## Function to make volcano plot
plot_volc <- function(top_genes, FDR_cut, model_name) {
    outGenes_plot <- top_genes %>%
        select(logFC, P.Value, adj.P.Val, ensemblID, Symbol)

    keyvals <- ifelse(
        outGenes_plot$P.Value > FDR_cut, "#f0e3d6", "#2a9d8f"
    )

    names(keyvals)[keyvals == "#2a9d8f"] <- paste0("pvalue < ", FDR_cut)
    names(keyvals)[keyvals == "#f0e3d6"] <- "Not significant"


    volcano_plot <- EnhancedVolcano(outGenes_plot,
        x = "logFC",
        y = "P.Value",
        selectLab = c(""),
        pCutoff = FDR_cut,
        FCcutoff = 0,
        lab = rownames(outGenes_plot),
        colCustom = keyvals
    ) + ylim(c(0, 5))

    pdf(paste0(out_plot, "/VolcanoPlot_", model_name, ".pdf"),
        height = 10,
        width = 8
    )
    print(volcano_plot)
    dev.off()
}

###############################################################################



##################################### DEA #####################################

## Model qc-snpPCs-Hb
## WITH -
## WITHOUT - totalAssignedGene, qSVs and tot.Thal
formula <- ~ PrimaryDx + AgeDeath + Flowcell + mitoRate + rRNA_rate + RIN + abs_ERCCsumLogErr +
    snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 +
    tot.Hb

res_formula <- DE_analysis(
    rse_gene = rse_gene,
    formula = formula,
    coef = "PrimaryDxSchizo",
    model_name = "qc-snpPCs-Hb"
)

plot_volc(res_formula, FDR_cut = 5e-02, model_name = "qc-snpPCs-Hb")

summary(res_formula$adj.P.Val)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.1267  0.1267  0.1267  0.1568  0.1408  0.9996


## Model qc-snpPCs-Hb-qSVs
## WITH - qSVs
## WITHOUT - totalAssignedGene and tot.Thal
formula <- ~ PrimaryDx + AgeDeath + Flowcell + mitoRate + rRNA_rate + RIN + abs_ERCCsumLogErr +
    snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 +
    qSV1 + qSV2 + qSV3 + qSV4 + qSV5 +
    tot.Hb

res_formula <- DE_analysis(
    rse_gene = rse_gene,
    formula = formula,
    coef = "PrimaryDxSchizo",
    model_name = "qc-snpPCs-Hb-qSVs"
)

plot_volc(
    res_formula,
    FDR_cut = 5e-02,
    model_name = "qc-snpPCs-Hb-qSVs"
)

summary(res_formula$adj.P.Val)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.3325  0.3325  0.3325  0.3704  0.3417  1.0000


## Model qc-totAssGene-snpPCs-Hb
## WITH - totalAssignedGene
## WITHOUT - qSVs and tot.Thal
formula <- ~ PrimaryDx + AgeDeath + Flowcell + mitoRate + rRNA_rate + RIN + totalAssignedGene + abs_ERCCsumLogErr +
    snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 +
    tot.Hb

res_formula <- DE_analysis(
    rse_gene = rse_gene,
    formula = formula, coef = "PrimaryDxSchizo",
    model_name = "qc-totAssGene-snpPCs-Hb"
)

plot_volc(res_formula,
    FDR_cut = 5e-02,
    model_name = "qc-totAssGene-snpPCs-Hb"
)

summary(res_formula$adj.P.Val)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.2338  0.2338  0.2338  0.2725  0.2541  0.9992


## Model qc-snpPCs-Hb-Thal
## WITH - tot.Thal
## WITHOUT - totalAssignedGene and qSVs
formula <- ~ PrimaryDx + AgeDeath + Flowcell + mitoRate + rRNA_rate + RIN + abs_ERCCsumLogErr +
    snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 +
    tot.Hb + tot.Thal

res_formula <- DE_analysis(
    rse_gene = rse_gene,
    formula = formula, coef = "PrimaryDxSchizo",
    model_name = "qc-snpPCs-Hb-Thal"
)

plot_volc(res_formula,
    FDR_cut = 5e-02,
    model_name = "qc-snpPCs-Hb-Thal"
)

summary(res_formula$adj.P.Val)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.2930  0.2930  0.2930  0.3270  0.3026  0.9982


## Model qc-totAssGene-snpPCs-Hb-qSVs
## WITH - qSVs and totalAssignedGene
## WITHOUT - tot.Thal
formula <- ~ PrimaryDx + AgeDeath + Flowcell + mitoRate + rRNA_rate + RIN + totalAssignedGene + abs_ERCCsumLogErr +
    snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 +
    qSV1 + qSV2 + qSV3 + qSV4 + qSV5 +
    tot.Hb

res_formula <- DE_analysis(
    rse_gene = rse_gene,
    formula = formula, coef = "PrimaryDxSchizo",
    model_name = "qc-totAssGene-snpPCs-Hb-qSVs"
)

plot_volc(res_formula,
    FDR_cut = 5e-02,
    model_name = "qc-totAssGene-snpPCs-Hb-qSVs"
)

summary(res_formula$adj.P.Val)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.4664  0.4664  0.4664  0.5008  0.4720  0.9993


## Model qc-snpPCs-Hb-Thal-qSVs
## WITH - qSVs and tot.Thal
## WITHOUT - totalAssignedGene
formula <- ~ PrimaryDx + AgeDeath + Flowcell + mitoRate + rRNA_rate + RIN + abs_ERCCsumLogErr +
    snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 +
    qSV1 + qSV2 + qSV3 + qSV4 + qSV5 +
    tot.Hb + tot.Thal

res_formula <- DE_analysis(
    rse_gene = rse_gene,
    formula = formula, coef = "PrimaryDxSchizo",
    model_name = "qc-snpPCs-Hb-Thal-qSVs"
)

plot_volc(res_formula,
    FDR_cut = 5e-02,
    model_name = "qc-snpPCs-Hb-Thal-qSVs"
)

summary(res_formula$adj.P.Val)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.5957  0.5957  0.5957  0.6236  0.5967  0.9999


## Model qc-totAssGene-snpPCs-Hb-Thal
## WITH - totalAssignedGene and tot.Thal
## WITHOUT - qSVs
formula <- ~ PrimaryDx + AgeDeath + Flowcell + mitoRate + rRNA_rate + RIN + totalAssignedGene + abs_ERCCsumLogErr +
    snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 +
    tot.Hb + tot.Thal

res_formula <- DE_analysis(
    rse_gene = rse_gene,
    formula = formula,
    coef = "PrimaryDxSchizo",
    model_name = "qc-totAssGene-snpPCs-Hb-Thal"
)

plot_volc(res_formula,
    FDR_cut = 5e-02,
    model_name = "qc-totAssGene-snpPCs-Hb-Thal"
)

summary(res_formula$adj.P.Val)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.3990  0.3990  0.3990  0.4345  0.4090  0.9997


## Model qc-totAssGene-snpPCs-Hb-Thal-qSVs
## WITH - totalAssignedGene, tot.Thal and qSVs
## WITHOUT -
formula <- ~ PrimaryDx + AgeDeath + Flowcell + mitoRate + rRNA_rate + RIN + totalAssignedGene + abs_ERCCsumLogErr +
    snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 +
    qSV1 + qSV2 + qSV3 + qSV4 + qSV5 +
    tot.Hb + tot.Thal

res_formula <- DE_analysis(
    rse_gene = rse_gene,
    formula = formula,
    coef = "PrimaryDxSchizo",
    model_name = "qc-totAssGene-snpPCs-Hb-Thal-qSVs"
)

plot_volc(res_formula,
    FDR_cut = 5e-02,
    model_name = "qc-totAssGene-snpPCs-Hb-Thal-qSVs"
)

summary(res_formula$adj.P.Val)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.6917  0.6917  0.6917  0.7142  0.6917  0.9998

###############################################################################



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
