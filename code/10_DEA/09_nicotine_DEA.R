#   For the second round of manuscript reviews, examine if nicotine toxicology
#   results confound case-control signal

library(here)
library(SummarizedExperiment)
library(edgeR)
library(limma)
library(tidyverse)
library(readxl)
library(qvalue)
library(sessioninfo)

rse_path = here('processed-data', 'rse_objects', 'rse_gene_Habenula_Pilot.rda')
supp_path = here('processed-data', '10_DEA', 'submission6_supp_tables.xlsx')
dx_dea_path = here(
    'processed-data', '10_DEA', '04_DEA', 
    'DEA_All-gene_qc-totAGene-qSVs-Hb-Thal.tsv'
)
plot_dir = here('plots', '10_DEA', '09_nicotine_DEA')
covars = c(
    'nicotine_tox', 'AgeDeath', 'Flowcell', 'mitoRate', 'rRNA_rate', 'RIN',
    'totalAssignedGene', 'abs_ERCCsumLogErr', 'qSV1', 'qSV2', 'qSV3', 'qSV4',
    'qSV5', 'qSV6', 'qSV7', 'qSV8', 'tot.Hb', 'tot.Thal'
)
coef_variable = 'nicotine_toxTRUE'
fdr = 0.1

dir.create(plot_dir, showWarnings = FALSE)

rse_gene = get(load(rse_path))
rse_gene = rse_gene[, rse_gene$PrimaryDx == 'Control']

supp_df = read_excel(supp_path, sheet = 'Table S1') |>
    select(BrNum, `Nicotine Tox`) |>
    rename(nicotine_tox = `Nicotine Tox`)

#   Merge in nicotine toxicology column to colData
col_data = colData(rse_gene) |>
    as_tibble() |>
    left_join(supp_df, by = 'BrNum') |>
    DataFrame()
rownames(col_data) = colnames(rse_gene)
colData(rse_gene) = col_data

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
    eb_feat, coef = coef_variable, xlab = "Mean of normalized counts",
    ylab = "logFC"
)
dev.off()

## Plot -log(p-value) vs logFC
pdf(file.path(plot_dir, 'volcano.pdf'))
volcanoplot(eb_feat, coef = coef_variable)
dev.off()

## Select top-ranked genes
top_genes = topTable(
        eb_feat, coef = coef_variable, p.value = 1, number = nrow(rse_scaled),
        sort.by = "none"
    ) |>
    as_tibble()

#   There are no DEGs
stopifnot(!any(top_genes$adj.P.Val < fdr))

## Histogram of p values
pdf(file.path(plot_dir, 'p_val_hist.pdf'))
hist(top_genes$P.Value, xlab = "p-value", main = "")
dev.off()

#   Read in case-control DEA results
dx_dea = read_tsv(dx_dea_path, show_col_types = FALSE)
stopifnot(identical(dx_dea$ensemblID, top_genes$ensemblID))

#   Compute replication statistics based on pi1. Among case-control DEGs, what
#   fraction are expected to be nicotine DEGs? Note there are no nicotine DEGs
#   to ask the reverse question

#   Compute pi1, adding a dummy p value of 1 if there's an error
p = top_genes$P.Value[dx_dea$adj.P.Val < fdr]
pi1 = tryCatch(
    1 - qvalue(p)$pi0, error = function(e) 1 - qvalue(c(1, p))$pi0
)

#   Gather nicotine and case-control t-stats
dea_df = tibble(
    gene = top_genes$ensemblID,
    t_nicotine = top_genes$t,
    t_dx = dx_dea$t
)

cor_obj = cor.test(dea_df$t_nicotine, dea_df$t_dx, method = "spearman")
anno_text = sprintf(
    " rho = %s\n p = %s\n Rep. y->x: %s%%\n",
    signif(cor_obj$estimate, 3),
    signif(cor_obj$p.value, 3),
    signif(100 * pi1, 3)
)

#   Is signal for nicotine status among controls correlated with
#   case-control signal?
p = ggplot(dea_df, aes(x = t_nicotine, y = t_dx)) +
    geom_bin2d(bins = 70) +
    geom_text(
        label = anno_text, x = -Inf, y = -Inf, hjust = 0, vjust = 0,
        size = 6
    ) +
    scale_fill_viridis_c() +
    labs(x = "Nicotine t-stat", y = "Case-control t-stat") +
    theme_bw(base_size = 15)
pdf(file.path(plot_dir, 'nicotine_vs_dx.pdf'))
print(p)
dev.off()

session_info()
