library(here)
library(tidyverse)
library(readxl)
library(SummarizedExperiment)
library(snpStats)

eqtl_path = here('processed-data', '17_eQTL', 'tensorQTL_output', 'FDR05.csv')
deg_path = here(
    'processed-data', '10_DEA', '04_DEA',
    'DEA_All-gene_qc-totAGene-qSVs-Hb-Thal.tsv'
)
gwas_path = here('processed-data', '17_eQTL', 'trubetskoy_gwas_supp_tab1.xls')
rse_path = here(
    'processed-data', 'rse_objects', 'rse_gene_Habenula_Pilot.rda'
)
plink_path_prefix = here(
    "processed-data", '08_bulk_snpPC', "habenula_genotypes"
)
plot_dir = here('plots', '17_eQTL')

sig_cutoff = 0.05
sig_cutoff_gwas = 5e-8

dir.create(plot_dir, showWarnings = FALSE)

#   Read in eQTL, DEA, and GWAS results as tibbles
eqtl = read_csv(eqtl_path, show_col_types = FALSE)
deg = read_tsv(deg_path, show_col_types = FALSE) |>
    filter(adj.P.Val < sig_cutoff)
gwas = read_excel(gwas_path) |>
    filter(P < sig_cutoff_gwas)

eqtl_gene = eqtl |>
    group_by(phenotype_id) |>
    summarize(n = n())

message(
    sprintf("%s significant SNP-gene pairs at FDR<%s", nrow(eqtl), sig_cutoff)
)
message(
    sprintf(
        "%s Genes total and median %s SNPs per gene",
        nrow(eqtl_gene),
        median(eqtl_gene$n)
    )
)

table(deg$gencodeID %in% eqtl_gene$phenotype_id)
deg$Symbol[deg$gencodeID %in% eqtl_gene$phenotype_id]

overlap_variants = eqtl |>
    filter(phenotype_id %in% deg$gencodeID) |>
    pull(variant_id)

rse_gene = get(load(rse_path))
colnames(rse_gene) = rse_gene$BrNum

#   Grab the expression for all genes significant in the DEA having significantly
#   associated variants found in the eQTL analysis; convert to long format
rse_small = rse_gene[
    deg$gencodeID[deg$gencodeID %in% eqtl_gene$phenotype_id],
]
express = assays(rse_small)$logcounts |>
    as.data.frame() |>
    rownames_to_column("gene_id") |>
    pivot_longer(cols = -gene_id, names_to = "sample_id", values_to = "logcount")

#   Read in genotype data and join with expression data
a = read.plink(paste0(plink_path_prefix, '.bed'))
exp_df = a$genotypes[, overlap_variants] |>
    as.data.frame() |>
    rownames_to_column("sample_id") |>
    pivot_longer(
        cols = -sample_id, names_to = "snp_id", values_to = "genotype"
    ) |>
    mutate(
        genotype = factor(as.integer(genotype)),
        gene_id = eqtl$phenotype_id[match(snp_id, eqtl$variant_id)]
    ) |>
    left_join(express, by = c('sample_id', 'gene_id')) |>
    filter(sample_id %in% express$sample_id)

#   Plot expression by genotype of each variant with one gene per page and
#   potentially several variants per gene
plot_list = list()
for (this_gene in unique(exp_df$gene_id)) {
    plot_list[[this_gene]] = exp_df |>
        filter(gene_id == this_gene) |>
        ggplot(mapping = aes(x = genotype, y = logcount)) +
            geom_boxplot(outlier.shape = NA) +
            geom_jitter() +
            facet_wrap(~ snp_id) +
            labs(
                title = rowData(rse_gene)$Symbol[
                    match(this_gene, rownames(rse_gene))
                ]
            )
}

pdf(file.path(plot_dir, 'expr_by_geno.pdf'))
print(plot_list)
dev.off()
