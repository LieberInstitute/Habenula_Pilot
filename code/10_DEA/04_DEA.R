library("here")
library("SummarizedExperiment")
library("edgeR")
library("limma")
library("dplyr")
library("gplots")
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


## Function to do all the DEA analysis with limma
DE_analysis <- function(rse_gene, formula, coef, model_name) {

    pdf(file = paste0(out_plot, "/DEAplots_", model_name, ".pdf"))
    # par(mfrow = c(3, 2))
    par(cex = 0.7, mai = c(0.1, 0.1, 0.1, 0.1))

    ## Model matrix
    model <- model.matrix(formula, data = colData(rse_gene))

    ## Use previous norm factors to scale the raw library sizes
    RSE_scaled <- calcNormFactors(rse_gene)

    par(fig = c(0.05, 0.5, 0.55, 0.95))
    ## Transform counts to log2(CPM): estimate mean-variance relationship for
    ## each gene
    vGene <- voom(RSE_scaled, design = model, plot = TRUE)

    ## Fit linear model for each gene
    fitGene <- lmFit(vGene)

    ## Empirical Bayesian calculation to obtain our significant genes: compute
    ## moderated F and t-statistics, and log-odds of DE
    eBGene <- eBayes(fitGene)

    par(fig = c(0.55, 1, 0.55, 0.95), new = TRUE)
    ## Plot average log expression vs logFC
    limma::plotMA(eBGene,
        coef = coef, xlab = "Mean of normalized counts",
        ylab = "logFC"
    )

    par(fig = c(0.05, 0.5, 0.1, 0.5), new = TRUE)

    ## Plot -log(p-value) vs logFC
    volcanoplot(eBGene, coef = coef)

    ## Select top-ranked genes
    top_genes <- topTable(eBGene, coef = coef, p.value = 1, number = nrow(rse_gene), sort.by = "none")

    ## Histogram of adjusted p values
    par(fig = c(0.55, 1, 0.1, 0.5), new = TRUE)
    hist(top_genes$P.Value, xlab = "p-value", main = "")

    par(fig = c(0.25, 0.75, 0, 0.1), new = TRUE)
    textplot(capture.output(summary(top_genes$adj.P.Val)))

    par(fig = c(0.0, 0.25, 0, 0.1), new = TRUE)
    textplot(capture.output(length(unique(top_genes$adj.P.Val))))

    dev.off()

    write.table(top_genes, file = paste0(out_data, "/DEA_AllGenes_", model_name, ".tsv"), sep = "\t", quote = FALSE)

    return(top_genes)
}

## Function to make volcano plot
plot_volc <- function(top_genes, FDR_cut, model_name) {
    outGenes_plot <- top_genes %>%
        select(logFC, P.Value, adj.P.Val, ensemblID, Symbol)

    keyvals <- ifelse(
        outGenes_plot$adj.P.Val > FDR_cut, "#f0e3d6", "#2a9d8f"
    )

    names(keyvals)[keyvals == "#2a9d8f"] <- paste0("FDR < ", FDR_cut)
    names(keyvals)[keyvals == "#f0e3d6"] <- "Not significant"


    volcano_plot <- EnhancedVolcano(outGenes_plot,
        x = "logFC",
        y = "adj.P.Val",
        selectLab = c(""),
        pCutoff = FDR_cut,
        FCcutoff = 0,
        lab = rownames(outGenes_plot),
        colCustom = keyvals
    ) + ylim(c(0, 2))

    pdf(paste0(out_plot, "/VolcanoPlot_", model_name, ".pdf"),
        height = 10,
        width = 8
    )
    print(volcano_plot)
    dev.off()
}

###############################################################################



##################################### DEA #####################################

## Model PrimaryDx
## WITH -
## WITHOUT - Everything
formula <- ~ PrimaryDx

res_formula <- DE_analysis(
    rse_gene = rse_gene,
    formula = formula,
    coef = "PrimaryDxSchizo",
    model_name = "PrimaryDx"
)

plot_volc(res_formula, FDR_cut = 10e-02, model_name = "PrimaryDx")

summary(res_formula$adj.P.Val)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.04559 0.72070 0.84753 0.81183 0.94172 0.99994


## Model PrimaryDx_Flow
## WITH - Flowcell
## WITHOUT - Almost Everything
formula <- ~ PrimaryDx + Flowcell

res_formula <- DE_analysis(
    rse_gene = rse_gene,
    formula = formula,
    coef = "PrimaryDxSchizo",
    model_name = "PrimaryDx_Flow"
)

plot_volc(res_formula, FDR_cut = 10e-02, model_name = "PrimaryDx_Flow")

summary(res_formula$adj.P.Val)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.05414 0.72771 0.85266 0.81553 0.94332 1.00000


## Model PrimaryDx_Flow_Age
## WITH - Flowcell and AgeDeath
## WITHOUT - Almost everything
formula <- ~ PrimaryDx + AgeDeath + Flowcell

res_formula <- DE_analysis(
    rse_gene = rse_gene,
    formula = formula,
    coef = "PrimaryDxSchizo",
    model_name = "PrimaryDx_Flow_Age"
)

plot_volc(res_formula, FDR_cut = 10e-02, model_name = "PrimaryDx_Flow_Age")

summary(res_formula$adj.P.Val)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.04104 0.73602 0.86147 0.81994 0.94316 0.99999


## Model qc
## WITH - All QCs
## WITHOUT - Almost everything
formula <- ~ PrimaryDx + AgeDeath + Flowcell + mitoRate + rRNA_rate + RIN + abs_ERCCsumLogErr

res_formula <- DE_analysis(
    rse_gene = rse_gene,
    formula = formula,
    coef = "PrimaryDxSchizo",
    model_name = "qc"
)

plot_volc(res_formula, FDR_cut = 10e-02, model_name = "qc")

summary(res_formula$adj.P.Val)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.01055 0.31640 0.56593 0.55644 0.80130 0.99998


## Model qc
## WITH - All QCs + 2 SNP PCs
## WITHOUT - Almost everything
formula <- ~ PrimaryDx + AgeDeath + Flowcell + mitoRate + rRNA_rate + RIN + abs_ERCCsumLogErr +
    snpPC1 + snpPC2

res_formula <- DE_analysis(
    rse_gene = rse_gene,
    formula = formula,
    coef = "PrimaryDxSchizo",
    model_name = "qc-2snpPCs"
)

plot_volc(res_formula, FDR_cut = 10e-02, model_name = "qc-2snpPCs")

summary(res_formula$adj.P.Val)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.01243 0.31284 0.56165 0.55427 0.79952 0.99998


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

plot_volc(res_formula, FDR_cut = 10e-02, model_name = "qc-snpPCs-Hb")

summary(res_formula$adj.P.Val)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.01578 0.32924 0.57105 0.57022 0.81831 0.99998


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
    FDR_cut = 10e-02,
    model_name = "qc-snpPCs-Hb-qSVs"
)

summary(res_formula$adj.P.Val)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.04399 0.45380 0.66999 0.65079 0.85941 0.99992


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
    FDR_cut = 10e-02,
    model_name = "qc-totAssGene-snpPCs-Hb"
)

summary(res_formula$adj.P.Val)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.009171 0.481927 0.697310 0.666823 0.870947 0.999976


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
    FDR_cut = 10e-02,
    model_name = "qc-snpPCs-Hb-Thal"
)

summary(res_formula$adj.P.Val)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.03448 0.57053 0.76840 0.72697 0.90941 0.99994


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
    FDR_cut = 10e-02,
    model_name = "qc-totAssGene-snpPCs-Hb-qSVs"
)

summary(res_formula$adj.P.Val)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.04834 0.56225 0.75461 0.72241 0.89503 0.99975


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
    FDR_cut = 10e-02,
    model_name = "qc-snpPCs-Hb-Thal-qSVs"
)

summary(res_formula$adj.P.Val)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.08682 0.76559 0.88306 0.84057 0.95633 1.00000



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
    FDR_cut = 10e-02,
    model_name = "qc-totAssGene-snpPCs-Hb-Thal"
)

summary(res_formula$adj.P.Val)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.01887 0.67658 0.83063 0.78162 0.92944 1.00000


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
    FDR_cut = 10e-02,
    model_name = "qc-totAssGene-snpPCs-Hb-Thal-qSVs"
)

summary(res_formula$adj.P.Val)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.07587 0.79389 0.90266 0.86291 0.96492 0.99996

###############################################################################



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
