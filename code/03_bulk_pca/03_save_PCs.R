library("SummarizedExperiment")
library("sessioninfo")
library("here")
library("sva")

rse_path = here(
    'processed-data', 'rse_objects', 'rse_gene_Habenula_Pilot.rda'
)
out_path = here('processed-data', '03_bulk_pca', 'PCs.rds')
covariates = c(
    "PrimaryDx", 'snpPC1', 'snpPC2', 'snpPC3', 'snpPC4', 'snpPC5'
)

rse = get(load(rse_path))

pd = as.data.frame(colData(rse)[, covariates])
mod <- model.matrix(
    as.formula(paste('~', paste(covariates, collapse = " + "))), data = pd
)

pca <- prcomp(t(assays(rse)$logcounts))
k <- num.sv(assays(rse)$logcounts, mod)

saveRDS(pca$x[, seq(k)], out_path)

session_info()
