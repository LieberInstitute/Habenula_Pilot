library(here)
library(tidyverse)
library(SummarizedExperiment)
library(sessioninfo)
library(data.table)
library(jaffelab)
library(cowplot)

eqtl_independent_path = here(
    'processed-data', '17_eQTL', 'tensorQTL_output', 'independent', 'FDR05.csv'
)
eqtl_int_path = here(
    'processed-data', '17_eQTL', 'tensorQTL_output', 'combined_interaction_subset.csv'
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
bsp2_dir = '/dcs05/lieber/liebercentral/BrainSEQ_LIBD001/brainseq_phase2/brainseq_phase2/browser'
bsp2_dlpfc_eqtl_path = file.path(bsp2_dir, 'BrainSeqPhaseII_eQTL_dlpfc_full.txt')
bsp2_hippo_eqtl_path = file.path(bsp2_dir, 'BrainSeqPhaseII_eQTL_hippo_full.txt')
bsp2_snp_path = file.path(bsp2_dir, 'BrainSeqPhaseII_snp_annotation.txt')
bsp2_geno_path = file.path(bsp2_dir, 'BrainSeqPhaseII_snp_genotype.txt')
bsp2_dlpfc_exp_path = file.path(bsp2_dir, 'BrainSeqPhaseII_clean_expression_eqtl_dlpfc_gene.txt')
bsp2_hippo_exp_path = file.path(bsp2_dir, 'BrainSeqPhaseII_clean_expression_eqtl_hippo_gene.txt')

plot_dir = here('plots', '17_eQTL')

sig_cutoff_deg = 0.1

################################################################################
#   Read in and preprocess data
################################################################################

#-------------------------------------------------------------------------------
#   DEGs and ~20k GWAS SNPs
#-------------------------------------------------------------------------------

#   FDR < 0.1 DEGs
deg = read_tsv(deg_path, show_col_types = FALSE) |>
    filter(adj.P.Val < sig_cutoff_deg)

#   ~20k "wide" GWAS SNPs
gwas_wide = read_csv(gwas_wide_path, show_col_types = FALSE)

#-------------------------------------------------------------------------------
#   Independent and interaction habenula eQTLs
#-------------------------------------------------------------------------------

#   Read in all significant independent eQTLs
eqtl_independent = read_csv(eqtl_independent_path, show_col_types = FALSE) |>
    mutate(pair_id = paste(phenotype_id, variant_id, sep = '_'))

#   The same eQTLs but taking stats from the interaction model
eqtl_int = read_csv(eqtl_int_path)

#   Warn that not all eQTLs were present in interaction models
message(
    sprintf(
        "Only %s of %s independent eQTLs observed in habenula interaction model",
        nrow(eqtl_int |> filter(interaction_var == "hb")),
        nrow(eqtl_independent)
    )
)
message(
    sprintf(
        "Only %s of %s independent eQTLs observed in thalamus interaction model",
        nrow(eqtl_int |> filter(interaction_var == "thal")),
        nrow(eqtl_independent)
    )
)

#-------------------------------------------------------------------------------
#   BSP2 eQTLs and genotypes
#-------------------------------------------------------------------------------

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
    filter(pair_id %in% eqtl_independent$pair_id)

bsp2_hippo = fread(bsp2_hippo_eqtl_path) |>
    as_tibble() |>
    filter(Type == 'gene') |>
    mutate(
        variant_id = bsp2_snp$variant_id[match(snp, bsp2_snp$snp)],
        pair_id = paste(feature_id, variant_id, sep = '_')
    ) |>
    filter(pair_id %in% eqtl_independent$pair_id)


#   Read in genotypes for SNPs that in habenula have eQTLs paired with habenula
#   DEGs
deg_variants = eqtl_independent |>
    filter(phenotype_id %in% deg$gencodeID) |>
    pull(variant_id)

bsp2_geno = fread(bsp2_geno_path) |>
    mutate(variant_id = bsp2_snp$variant_id[match(V1, bsp2_snp$snp)]) |>
    filter(variant_id %in% deg_variants)

bsp2_exp_dlpfc = fread(bsp2_dlpfc_exp_path)

#   Warn that not all eQTLs were present in BSP2
message(
    sprintf(
        "Only %s of %s independent eQTLs observed in BSP2 DLPFC",
        nrow(bsp2_dlpfc),
        nrow(eqtl_independent)
    )
)
message(
    sprintf(
        "Only %s of %s independent eQTLs observed in BSP2 hippo",
        nrow(bsp2_hippo),
        nrow(eqtl_independent)
    )
)

################################################################################
#   Thalamus vs. habenula interaction plots for eQTLs overlapping DEGs or
#   wide GWAS SNPs
################################################################################

#   Keep a filtered copy to those overlapping either the wide GWAS or FDR < 0.1
#   DEGs
filt_eqtl_independent = eqtl_independent |>
    filter(
        phenotype_id %in% deg$gencodeID | 
        variant_id %in% gwas_wide$variant_id
    )

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

hb_pairs = eqtl_int |> filter(interaction_var == "hb") |> pull(pair_id)
thal_pairs = eqtl_int |> filter(interaction_var == "thal") |> pull(pair_id)
eqtl_int_both = eqtl_int |>
    dplyr::filter(
        pair_id %in% hb_pairs,
        pair_id %in% thal_pairs,
        pair_id %in% filt_eqtl_independent$pair_id
    )

#   Make note of filtered independent pairs not measured in each interaction
#   model
message(
    sprintf(
        "Pairs present in independent but absent from at least one interaction model: %s",
        paste(
            setdiff(filt_eqtl_independent$pair_id, eqtl_int_both$pair_id),
            collapse = ', '
        )
    )
)

p = eqtl_int_both |>
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
        geom_point(size = 3) +
        geom_hline(yintercept = 0) +
        geom_vline(xintercept = 0) +
        theme_bw(base_size = 20) +
        labs(
            x = "Beta: Habenula Interaction", y = "Beta: Thalamus Interaction",
            color = "eQTL Overlap"
        )

pdf(file.path(plot_dir, 'interaction_beta.pdf'), width = 10, height = 7)
print(p)
dev.off()

################################################################################
#   BSP2 vs habenula beta value plots at the same eQTLs
################################################################################

p_dlpfc = bsp2_dlpfc |>
    select(pair_id, beta) |>
    dplyr::rename(beta_DLPFC = beta) |>
    left_join(
        eqtl_independent |>
            select(pair_id, slope) |>
            dplyr::rename(beta_habenula = slope),
        by = "pair_id"
    ) |>
    ggplot(mapping = aes(x = beta_habenula, y = beta_DLPFC)) +
        geom_point() +
        geom_abline(slope = 1) +
        theme_bw(base_size = 20) +
        labs(x = "Beta: Habenula", y = "Beta: DLPFC")

p_hippo = bsp2_hippo |>
    select(pair_id, beta) |>
    dplyr::rename(beta_hippo = beta) |>
    left_join(
        eqtl_independent |>
            select(pair_id, slope) |>
            dplyr::rename(beta_habenula = slope),
        by = "pair_id"
    ) |>
    ggplot(mapping = aes(x = beta_habenula, y = beta_hippo)) +
        geom_point() +
        geom_abline(slope = 1) +
        theme_bw(base_size = 20) +
        labs(x = "Beta: Habenula", y = "Beta: Hippocampus")

pdf(file.path(plot_dir, 'BSP2_vs_habenula_beta.pdf'), width = 12, height = 6)
print(plot_grid(plotlist = list(p_dlpfc, p_hippo)))
dev.off()

pdf(
    file.path(plot_dir, 'BSP2_vs_habenula_beta_coord_fixed.pdf'),
    width = 12, height = 6
)
print(
    plot_grid(plotlist = list(p_dlpfc + coord_fixed(), p_hippo + coord_fixed()))
)
dev.off()

################################################################################
#   BSP2 genotypes at DEG-paired significant independent habenula eQTLs
################################################################################

session_info()
