library(here)
library(tidyverse)
library(readxl)

eqtl_path = here('processed-data', '17_eQTL', 'tensorQTL_output', 'FDR05.csv')
deg_path = here(
    'processed-data', '10_DEA', '04_DEA',
    'DEA_All-gene_qc-totAGene-qSVs-Hb-Thal.tsv'
)
gwas_path = here('processed-data', '17_eQTL', 'trubetskoy_gwas_supp_tab1.xls')

sig_cutoff = 0.05
sig_cutoff_gwas = 5e-8

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
