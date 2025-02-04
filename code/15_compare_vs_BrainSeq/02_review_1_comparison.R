#   After the first round of manuscript reviews, we were asked to perform
#   threshold-free comparisons of DE results between habenula and other
#   brain regions. This script performs these comparisons

library(here)
library(tidyverse)
library(BiocFileCache)
library(readxl)
library(sessioninfo)
library(cowplot)

plot_dir = here("plots", "15_compare_vs_BrainSeq", "review_round_1")
hab_de_path = here(
    "processed-data", "10_DEA", "04_DEA",
    'DEA_All-gene_qc-totAGene-qSVs-Hb-Thal.tsv'
)
bsp2_dlpfc_de_path = here(
    "raw-data", "15_compare_vs_BrainSeq",
    "dxStats_dlpfc_filtered_qSVA_geneLevel_noHGoldQSV_matchDLPFC.rda"
)
bsp2_hippo_de_path = here(
    "raw-data", "15_compare_vs_BrainSeq",
    "dxStats_hippo_filtered_qSVA_geneLevel_noHGoldQSV_matchHIPPO.rda"
)
bsp3_caudate_de_url = "https://caudate-eqtl.s3.us-west-2.amazonaws.com/BrainSeq_Phase3_Caudate_DifferentialExpression_DxSZ_all.txt.gz"
dg_de_path = here(
    "raw-data", "15_compare_vs_BrainSeq", "41593_2020_604_MOESM2_ESM.xlsx"
)
other_brain_regions = c("DLPFC", "Hippo", "Caudate", "DG")

dir.create(plot_dir, showWarnings = FALSE)

################################################################################
#   Functions
################################################################################

#   Given two tibbles of DE results, return a tibble containing
#   concordance-at-the-top results
to_cat_df = function(hab_de, other_de, region, max_x, num_samples) {
    x_vals = as.integer(round(seq_len(num_samples) * max_x / num_samples))

    y_vals = sapply(
        x_vals,
        function(x) {
            top_hab = hab_de |>
                arrange(adj.P.Val) |>
                slice_head(n = x) |>
                pull(ensemblID)
            top_other = other_de |> 
                arrange(adj.P.Val) |>
                slice_head(n = x) |>
                pull(ensemblID)

            return(length(intersect(top_hab, top_other)) / x)
        }
    )

    cat_df = tibble(x = x_vals, y = y_vals, region = region)
    return(cat_df)
}

#   Given two tibbles of DE results, return a scatterplot of t-statistics
#   across brain regions (habenula and 'region')
t_stat_plot = function(hab_de, other_de, region) {
    shared_de = inner_join(hab_de, other_de, by = "ensemblID")

    cor_obj = cor.test(shared_de$t.x, shared_de$t.y, method = "spearman")
    anno_text = sprintf(
        "rho = %s \np = %s \n",
        signif(cor_obj$estimate, 3),
        signif(cor_obj$p.value, 3)
    )
    
    p = ggplot(shared_de, aes(x = t.x, y = t.y)) +
        geom_bin2d(bins = 70) +
        geom_text(
            label = anno_text, x = Inf, y = -Inf, hjust = 1, vjust = 0,
            size = 4
        ) +
        scale_fill_viridis_c() +
        labs(x = "Habenula t-stat", y = paste(region, "t-stat")) +
        theme_bw(base_size = 15) +
        theme(legend.key.size = unit(0.35, 'cm'))

    return(p)
}

################################################################################
#   Load and clean DE results for all brain regions
################################################################################

de_list = list()

#   Habenula
de_list[['habenula']] = read_tsv(hab_de_path, show_col_types = FALSE) |>
    select(ensemblID, t, P.Value, adj.P.Val)

#   BSP2 DLPFC
load(bsp2_dlpfc_de_path, verbose = TRUE)
de_list[['DLPFC']] = outGene |>
    as_tibble() |>
    select(ensemblID, t, P.Value, adj.P.Val)
rm(outGene0, outGeneNoAdj, outGene)

#   BSP2 Hippocampus
load(bsp2_hippo_de_path, verbose = TRUE)
de_list[['Hippo']] = outGene |>
    as_tibble() |>
    select(ensemblID, t, P.Value, adj.P.Val)
rm(outGene0, outGeneNoAdj, outGene)

#   BSP3 Caudate
bfc <- BiocFileCache()
bsp3_caudate_file <- bfcrpath(
    bfc, bsp3_caudate_de_url, exact = TRUE
)
de_list[['Caudate']] = read_tsv(bsp3_caudate_file, show_col_types = FALSE) |>
    filter(Type == "Gene") |>
    select(ensemblID, t, P.Value, adj.P.Val)

#   DG
de_list[['DG']] = read_xlsx(dg_de_path, sheet = "Table_S10") |>
    select(ensemblID, SZ_t, SZ_P.Value, SZ_adj.P.Val) |>
    rename(t = SZ_t, P.Value = SZ_P.Value, adj.P.Val = SZ_adj.P.Val)

################################################################################
#   Concordance-at-the-top plots
################################################################################

cat_df_list = list()
for (region in other_brain_regions) {
    cat_df_list[[region]] = to_cat_df(
        hab_de = de_list[['habenula']],
        other_de = de_list[[region]],
        region = region,
        max_x = 3000,
        num_samples = 100
    )
}

p = do.call(rbind, cat_df_list) |>
    ggplot(aes(x = x, y = y)) +
    geom_line() +
    facet_wrap(~region) +
    labs(
        x = "Number of Top DE Genes",
        y = "Concordance",
        color = "Brain region"
    ) +
    theme_bw(base_size = 20)

pdf(file.path(plot_dir, "CAT_plots.pdf"))
print(p)
dev.off()

################################################################################
#   Concordance of t-statistics
################################################################################

plot_list = list()
for (region in other_brain_regions) {
    plot_list[[region]] = t_stat_plot(
        hab_de = de_list[['habenula']],
        other_de = de_list[[region]],
        region = region
    )
}
pdf(file.path(plot_dir, "t_stat_comparison.pdf"), width = 10)
plot_grid(plotlist = plot_list, ncol = 2)
dev.off()

session_info()
