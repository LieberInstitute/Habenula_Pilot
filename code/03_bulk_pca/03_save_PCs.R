library("SummarizedExperiment")
library("sessioninfo")
library("here")
library("sva")
library("tidyverse")

rse_path = here(
    'processed-data', 'rse_objects', 'rse_gene_Habenula_Pilot.rda'
)
out_path = here('processed-data', '03_bulk_pca', 'PCs.rds')
covariates = c(
    "PrimaryDx", 'snpPC1', 'snpPC2', 'snpPC3', 'snpPC4', 'snpPC5'
)
snp_pcs_path = here(
    'processed-data', '08_bulk_snpPC',
    'Hb_gt_merged_R.9_MAF.05_ann_filt.snpPCs.tab'
)

set.seed(0)

#   Load RSE and overwrite SNP PCs with the most recent values computed from the
#   properly filtered genotyping data
rse = get(load(rse_path))
colData(rse) = colData(rse) |>
    as_tibble() |>
    select(!matches('^snpPC')) |>
    left_join(
        read_tsv(snp_pcs_path, show_col_types = FALSE) |>
            #   Fix the name of one donor
            mutate(BrNum = ifelse(BrNum == "Br0983", "Br983", BrNum)),
        by = 'BrNum'
    ) |>
    DataFrame()

pd = as.data.frame(colData(rse)[, covariates])
mod <- model.matrix(
    as.formula(paste('~', paste(covariates, collapse = " + "))), data = pd
)

pca <- prcomp(t(assays(rse)$logcounts))
k <- num.sv(assays(rse)$logcounts, mod)

pc_df = pca$x[, seq(k)] |>
    as_tibble() |>
    mutate(RNum = rse$RNum)
saveRDS(pc_df, out_path)

session_info()
