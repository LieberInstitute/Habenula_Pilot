library(here)
library(tidyverse)
library(SummarizedExperiment)
library(WGCNA)
library(jaffelab)

rse_path = here(
    'processed-data', 'rse_objects', 'rse_gene_Habenula_Pilot.rda'
)
deg_covariates = c(
    'PrimaryDx', 'AgeDeath', 'Flowcell', 'mitoRate', 'rRNA_rate', 'RIN',
    'totalAssignedGene', 'abs_ERCCsumLogErr', 'qSV1', 'qSV2', 'qSV3', 'qSV4',
    'qSV5', 'qSV6', 'qSV7', 'qSV8', 'tot.Hb', 'tot.Thal'
)

rse_gene = get(load(rse_path))

mod_deg <- model.matrix(
    as.formula(paste('~', paste(deg_covariates, collapse = " + "))),
    data = colData(rse_gene)
)

exp_mat = t(cleaningY(assays(rse_gene)$logcounts, mod_deg, P = 3))
