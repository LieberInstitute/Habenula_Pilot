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

bsp2_dir = '/dcs05/lieber/liebercentral/BrainSEQ_LIBD001/brainseq_phase2/brainseq_phase2/browser'
bsp2_dlpfc_eqtl_path = file.path(bsp2_dir, 'BrainSeqPhaseII_eQTL_dlpfc_full.txt')
bsp2_hippo_eqtl_path = file.path(bsp2_dir, 'BrainSeqPhaseII_eQTL_hippo_full.txt')
bsp2_snp_path = file.path(bsp2_dir, 'BrainSeqPhaseII_snp_annotation.txt')

bsp2_out_path = here('processed-data', '17_eQTL', 'BSP2_cleaned_expr_geno.csv')

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

#-------------------------------------------------------------------------------
#   BSP2 eQTLs with genotypes
#-------------------------------------------------------------------------------

#   Table of alternative SNP names
bsp2_snp = fread(bsp2_snp_path) |>
    as_tibble() |>
    #   Use the type of SNP ID used in the habenula analysis
    mutate(
        variant_id = sprintf(
            '%s:%s:%s:%s', chr_hg38, pos_hg38, newref, newcount
        )
    )

bsp2_dlpfc = fread(bsp2_dlpfc_eqtl_path) |>
    as_tibble() |>
    filter(Type == 'gene') |>
    mutate(
        variant_id = bsp2_snp$variant_id[match(snp, bsp2_snp$snp)],
        pair_id = paste(feature_id, variant_id, sep = '_')
    ) |>
    filter(pair_id %in% eqtl_independent$pair_id) |>
    mutate(region = "DLPFC")

bsp2_hippo = fread(bsp2_hippo_eqtl_path) |>
    as_tibble() |>
    filter(Type == 'gene') |>
    mutate(
        variant_id = bsp2_snp$variant_id[match(snp, bsp2_snp$snp)],
        pair_id = paste(feature_id, variant_id, sep = '_')
    ) |>
    filter(pair_id %in% eqtl_independent$pair_id) |>
    mutate(region = "Hippocampus")

bsp2 = rbind(bsp2_dlpfc, bsp2_hippo) |>
    dplyr::rename(gene_id = feature_id) |>
    select(gene_id, variant_id, pair_id, FDR, beta, region)

write_csv(bsp2, bsp2_out_path)

session_info()
