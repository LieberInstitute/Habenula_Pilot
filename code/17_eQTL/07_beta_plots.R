library(here)
library(tidyverse)
library(SummarizedExperiment)
library(sessioninfo)
library(data.table)
library(jaffelab)

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

#   Cache some temporary results from this script for quicker interactive testing
eqtl_int_path = here(
    'processed-data', '17_eQTL', 'tensorQTL_output', 'combined_interaaction_subset.csv'
)

deg_path = here(
    'processed-data', '10_DEA', '04_DEA',
    'DEA_All-gene_qc-totAGene-qSVs-Hb-Thal.tsv'
)
gwas_wide_path = here(
    'processed-data', '17_eQTL', 'gwas_wide_filtered.csv'
)
rse_path = here(
    'processed-data', 'rse_objects', 'rse_gene_Habenula_Pilot.rda'
)

plot_dir = here('plots', '17_eQTL')

sig_cutoff_deg = 0.1

deg = read_tsv(deg_path, show_col_types = FALSE) |>
    filter(adj.P.Val < sig_cutoff_deg)

gwas_wide = read_csv(gwas_wide_path, show_col_types = FALSE)

#   Read in all significant independent eQTLs and keep a filtered copy to those
#   overlapping either the wide GWAS or FDR < 0.1 DEGs
eqtl_independent = read_csv(eqtl_independent_path, show_col_types = FALSE)
filt_eqtl_independent = eqtl_independent |>
    filter(
        phenotype_id %in% deg$gencodeID | 
        variant_id %in% gwas_wide$variant_id
    ) |>
    mutate(pair_id = paste(phenotype_id, variant_id, sep = '_'))

#   Note that there are 11 SNPs overlapping the GWAS data, but these SNPs appear
#   in 13 eQTLs total
message(
    "Num SNPs overlapping wide GWAS data: ",
    eqtl_independent |>
        filter(variant_id %in% gwas_wide$variant_id) |>
        pull(variant_id) |>
        unique() |>
        length()
)
message(
    "Num eQTLs including those SNPs: ",
    eqtl_independent |>
        filter(variant_id %in% gwas_wide$variant_id) |>
        nrow()
)

#   Read in the full set of interaction eQTLs (with habenula and thalamus
#   fraction) and individually filter SNPs that overlap in the independent set
#   with the GWAS or DEGs (filtered one at a time to reduce memory usage)
eqtl_hb = read_csv(eqtl_hb_path, show_col_types = FALSE) |>
    mutate(
        interaction_var = "hb",
        pair_id = paste(phenotype_id, variant_id, sep = '_')
    ) |>
    filter(pair_id %in% filt_eqtl_independent$pair_id)
eqtl_thal = read_csv(eqtl_thal_path, show_col_types = FALSE) |>
    mutate(
        interaction_var = "thal",
        pair_id = paste(phenotype_id, variant_id, sep = '_')
    ) |>
    filter(pair_id %in% filt_eqtl_independent$pair_id)

eqtl_int = rbind(eqtl_hb, eqtl_thal) |>
    filter(pair_id %in% intersect(eqtl_hb$pair_id, eqtl_thal$pair_id))

write_csv(eqtl_int, eqtl_int_path)

#   For quicker interactive testing
# eqtl_int = read_csv(eqtl_int_path)

#   Make note of filtered independent pairs not measured in each interaction
#   model
message(
    sprintf(
        "Pairs present in independent but absent from at least one interaction model: %s",
        paste(
            setdiff(filt_eqtl_independent$pair_id, eqtl_int$pair_id),
            collapse = ', '
        )
    )
)

p = eqtl_int |>
    select(phenotype_id, variant_id, b_gi, pval_gi, interaction_var) |>
    mutate(
        source = factor(
            ifelse(variant_id %in% gwas_wide$variant_id, "GWAS SNP", "DEG")
        )
    ) |>
    pivot_wider(
        values_from = c("b_gi", "pval_gi"), names_from = "interaction_var"
    ) |>
    ggplot(mapping = aes(x = b_gi_hb, y = b_gi_thal, color = source)) +
        geom_point() +
        geom_hline(yintercept = 0) +
        geom_vline(xintercept = 0)

pdf(file.path(plot_dir, 'interaction_beta.pdf'))
print(p)
dev.off()
