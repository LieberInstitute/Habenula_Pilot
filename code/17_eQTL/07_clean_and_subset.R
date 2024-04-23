#   Many of the plots for 08_beta_plots.* require loading very large datasets,
#   only to use a tiny subset of the data. This script contains this
#   preprocessing such that the plotting script can simply load the tiny subset
#   and plot

library(here)
library(tidyverse)
library(SummarizedExperiment)
library(sessioninfo)
library(data.table)

eqtl_independent_path = here(
    'processed-data', '17_eQTL', 'tensorQTL_output', 'independent', 'FDR05.csv'
)
eqtl_hb_path = here(
    'processed-data', '17_eQTL', 'tensorQTL_output', 'interaction_tot_Hb',
    'all.csv'
)
eqtl_thal_path = here(
    'processed-data', '17_eQTL', 'tensorQTL_output', 'interaction_tot_Thal',
    'all.csv'
)
eqtl_int_path = here(
    'processed-data', '17_eQTL', 'tensorQTL_output', 'combined_interaction_subset.csv'
)

#-------------------------------------------------------------------------------
#   Interaction habenula eQTLs
#-------------------------------------------------------------------------------

#   Read in all significant independent eQTLs
eqtl_independent = read_csv(eqtl_independent_path, show_col_types = FALSE) |>
    mutate(pair_id = paste(phenotype_id, variant_id, sep = '_'))

#   Read in the full set of interaction eQTLs (with habenula and thalamus
#   fraction) and filter to those significant in independent model
eqtl_hb = fread(eqtl_hb_path) |>
    as_tibble() |>
    mutate(
        interaction_var = "hb",
        pair_id = paste(phenotype_id, variant_id, sep = '_')
    ) |>
    filter(pair_id %in% eqtl_independent$pair_id)

eqtl_thal = fread(eqtl_thal_path) |>
    as_tibble() |>
    mutate(
        interaction_var = "thal",
        pair_id = paste(phenotype_id, variant_id, sep = '_')
    ) |>
    filter(pair_id %in% eqtl_independent$pair_id)

rbind(eqtl_hb, eqtl_thal) |> write_csv(eqtl_int_path)
