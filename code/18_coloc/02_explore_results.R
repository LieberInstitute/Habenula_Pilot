library(here)
library(tidyverse)
library(sessioninfo)
library(SummarizedExperiment)
library(MRutils)

results_path = here("processed-data", "18_coloc", "coloc_results.rds")
out_path = here("processed-data", "18_coloc", "filtered_results.csv")
deg_path = here(
    'processed-data', '10_DEA', '04_DEA',
    'DEA_All-gene_qc-totAGene-qSVs-Hb-Thal.tsv'
)
eqtl_path = here(
    'processed-data', '17_eQTL', 'tensorQTL_output', 'independent', 'FDR05.csv'
)
rse_path = here(
    'processed-data', 'rse_objects', 'rse_gene_Habenula_Pilot.rda'
)
eqtl_overlap_path = here(
    "processed-data", "18_coloc", "eqtl_overlap.csv"
)
gwas_wide_filt_path = here(
    'processed-data', '17_eQTL', 'gwas_wide_filtered.csv'
)
supp_tab_path = here(
    "processed-data", "18_coloc", "supp_table.csv"
)

h4_cutoff = 0.8
deg_sig_cutoff = 0.1

results = readRDS(results_path)

shared_causal_variant = sapply(
    results, function(x) x$summary[['PP.H4.abf']] > h4_cutoff
)
message(sprintf("Genes with a shared causal variant (cutoff > %s):", h4_cutoff))
table(shared_causal_variant)

all_genes = names(results)

#   Add gene to each tibble of results
results_df = lapply(
    1:length(results),
    function(i) {
        as_tibble(results[[i]]$results) |>
            mutate(gene = all_genes[i])
    }
)

results_df = do.call(rbind, results_df) |>
    #   Only take genes with sufficient probability of shared causal variant
    filter(gene %in% all_genes[shared_causal_variant]) |>
    #   For each gene, take just the set of SNPs in the 95% credible set
    group_by(gene) |>
    arrange(desc(SNP.PP.H4)) |>
    mutate(h4_cumsum = cumsum(SNP.PP.H4)) |>
    filter(h4_cumsum <= 0.95) |>
    ungroup() |>
    #   Save just the gene, SNP, and shared-causal-variant probability
    select(gene, snp, SNP.PP.H4)

write_csv(results_df, out_path)

#   Load RSE just to get gene symbols
rse_gene = get(load(rse_path))

#   Print symbols of genes that colocalize
coloc_genes = rowData(rse_gene)$Symbol[
    match(unique(results_df$gene), rownames(rse_gene))
]
message(
    'Names of colocalized genes: "',
    paste(coloc_genes, collapse = '", "'),
    '"'
)

#   Check if genes overlap with DEGs
deg = read_tsv(deg_path, show_col_types = FALSE) |>
    filter(adj.P.Val < deg_sig_cutoff)

results_df$is_gene_de = results_df$gene %in% deg$gencodeID
message(
    sprintf(
        "Any coloc genes differentially expressed (FDR < %s)?", deg_sig_cutoff
    )
)
table(results_df$is_gene_de)

#   Check for overlaps with eQTLs
eqtl = read_csv(eqtl_path, show_col_types = FALSE)

results_df$is_eqtl = paste(results_df$gene, results_df$snp, sep = '_') %in%
    paste(eqtl$phenotype_id, eqtl$variant_id, sep = '_')

message("Any coloc gene-SNP pairs also independent eQTLs?")
table(results_df$is_eqtl)

#   Finally check SNP overlap with PGC3 risk SNPs
gwas = read_csv(gwas_wide_filt_path, show_col_types = FALSE)

results_df$is_risk_snp = results_df$snp %in% gwas$variant_id
message("Are any coloc SNPs PGC3 risk SNPs?")
table(results_df$is_risk_snp)

write_csv(results_df |> select(-SNP.PP.H4), supp_tab_path)

#   Write gene and SNP info for the colocalized gene-SNP pairs that overlap with
#   independent eQTLs
eqtl |>
    filter(
        paste(phenotype_id, variant_id, sep = '_') %in%
        paste(results_df$gene, results_df$snp, sep = '_')
    ) |>
    mutate(
        snp_rs_id = sapply(
            variant_id,
            function(x) {
                get_rsid_from_position(
                    chrom = str_split_i(x, ':', 1),
                    pos = as.integer(str_split_i(x, ':', 2)),
                    ref = str_split_i(x, ':', 3),
                    alt = str_split_i(x, ':', 4),
                    assembly = "hg38"
                )
            }
        ),
        gene_symbol = rowData(rse_gene)$Symbol[
            match(phenotype_id, rownames(rse_gene))
        ]
    ) |>
    dplyr::rename(gene_id = phenotype_id, snp_name = variant_id) |>
    dplyr::select(gene_id, gene_symbol, snp_name, snp_rs_id) |>
    write_csv(eqtl_overlap_path)

session_info()
