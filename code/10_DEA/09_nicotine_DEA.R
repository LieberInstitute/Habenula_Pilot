library(here)
library(SummarizedExperiment)
library(edgeR)
library(limma)
library(tidyverse)
library(sessioninfo)

rse_path = here('processed-data', 'rse_objects', 'rse_gene_Habenula_Pilot.rda')
plot_dir = here('plots', '10_DEA', '09_nicotine_DEA')
covars = c(
    'AgeDeath', 'Flowcell', 'mitoRate', 'rRNA_rate', 'RIN', 'totalAssignedGene',
    'abs_ERCCsumLogErr', 'qSV1', 'qSV2', 'qSV3', 'qSV4', 'qSV5', 'qSV6', 'qSV7',
    'qSV8', 'tot.Hb', 'tot.Thal'
)

dir.create(plot_dir, showWarnings = FALSE)

rse_gene = get(load(rse_path))
stopifnot(all(covars %in% colnames(colData(rse_gene))))

model = paste('~', paste(covars, collapse = " + ")) |>
    as.formula() |>
    model.matrix(data = colData(rse_gene))

## Use previous norm factors to scale the raw library sizes
rse_scaled = calcNormFactors(rse_gene)

## Transform counts to log2(CPM): estimate mean-variance relationship for
## each gene
pdf(file.path(plot_dir, 'voom.pdf'))
v_feat = voom(rse_scaled, design = model, plot = TRUE)
dev.off()

## Fit linear model for each gene
fit_feat = lmFit(v_feat)

## Empirical Bayesian calculation to obtain our significant genes: compute
## moderated F and t-statistics, and log-odds of DE
eb_feat = eBayes(fit_feat)

## Plot average log expression vs logFC
pdf(file.path(plot_dir, 'MA.pdf'))
plotMA(
    eb_feat, coef = covars[1], xlab = "Mean of normalized counts",
    ylab = "logFC"
)
dev.off()

## Plot -log(p-value) vs logFC
pdf(file.path(plot_dir, 'volcano.pdf'))
volcanoplot(eb_feat, coef = covars[1])
dev.off()

## Select top-ranked genes
top_genes = topTable(
        eb_feat, coef = covars[1], p.value = 1, number = nrow(rse_gene),
        sort.by = "none"
    ) |>
    as_tibble()

## Histogram of p values
pdf(file.path(plot_dir, 'p_val_hist.pdf'))
hist(top_genes$P.Value, xlab = "p-value", main = "")
dev.off()

session_info()
