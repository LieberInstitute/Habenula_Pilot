library(here)
library(tidyverse)
library(readxl)
library(SummarizedExperiment)
library(snpStats)
library(bigsnpr)
library(sessioninfo)
library(cowplot)

eqtl_path = here(
    'processed-data', '17_eQTL', 'tensorQTL_output', 'nominal', 'FDR05.csv'
)
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
paired_variants_path = here(
    "processed-data", "17_eQTL", "DEA_paired_variants.txt"
)
plot_dir = here('plots', '17_eQTL')

sig_cutoff = 0.05
sig_cutoff_gwas = 5e-8

covariates = c(
    "PrimaryDx", 'snpPC1', 'snpPC2', 'snpPC3', 'snpPC4', 'snpPC5'
)

lift_over_path = system('which liftOver', intern = TRUE)
dir.create(plot_dir, showWarnings = FALSE)

################################################################################
#   Read in eQTL, DEA, and GWAS results, and the RSE
################################################################################

eqtl = read_csv(eqtl_path, show_col_types = FALSE)

deg = read_tsv(deg_path, show_col_types = FALSE) |>
    filter(adj.P.Val < sig_cutoff)

rse_gene = get(load(rse_path))
colnames(rse_gene) = rse_gene$BrNum

gwas = read_excel(gwas_path) |>
    filter(P < sig_cutoff_gwas) |>
    dplyr::rename(chr = CHR, pos = BP) |>
    #   Map from hg19 to hg38 and drop anything that fails
    snp_modifyBuild(lift_over_path, from = 'hg19', to = 'hg38') |>
    filter(!is.na(pos)) |>
    #   Construct variant_id from SNP info
    mutate(
        variant_id = sprintf(
            'chr%s:%s:%s:%s',
            chr,
            pos,
            str_split_i(A1A2, '/', 1),
            str_split_i(A1A2, '/', 2)
        )
    )

################################################################################
#   Compare each type of result
################################################################################

#-------------------------------------------------------------------------------
#   Info about eQTLS alone
#-------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------
#   Compare eQTLs with GWAS SNPs
#-------------------------------------------------------------------------------

message("Number of Trubetskoy et al. GWAS SNPs matching sig eQTLs:")
table(gwas$variant_id %in% eqtl$variant_id)

gwas_paired_genes = eqtl |>
    filter(variant_id %in% gwas$variant_id) |>
    pull(phenotype_id)

message("Genes paired with GWAS-significant eQTL SNPs:")
rowData(rse_gene)$Symbol[
    match(gwas_paired_genes, rownames(rse_gene))
] |>
    paste(collapse = ', ') |>
    message()

#-------------------------------------------------------------------------------
#   Compare eQTLS with DE genes
#-------------------------------------------------------------------------------

message("Number and name of DE genes implicated in any eQTLs:")
table(deg$gencodeID %in% eqtl_gene$phenotype_id)
deg$Symbol[deg$gencodeID %in% eqtl_gene$phenotype_id] |>
    paste(collapse = ', ') |>
    message()

message("Genes implicated in the GWAS, DEA, and eQTLs:")
dea_paired_genes = deg$gencodeID[deg$gencodeID %in% eqtl_gene$phenotype_id]
dea_paired_genes[dea_paired_genes %in% gwas_paired_genes] |>
    paste(collapse = ', ') |>
    message()

#   Write paired variants to a text file. This will be read used to subset the
#   big VCF to ensure the below method for reading in genotypes (reading in the
#   plink bed file) works as VariantAnnotation::readVcf() does
dea_paired_variants = eqtl |>
    filter(phenotype_id %in% deg$gencodeID) |>
    pull(variant_id)
writeLines(dea_paired_variants, paired_variants_path)

################################################################################
#   Plots exploring how genotype affects expression at select eQTLs
################################################################################

#   Model covariates to later call jaffelab::cleaningY() to residualize
#   expression for downstream plots
pd = as.data.frame(colData(rse_gene)[, covariates])
mod <- model.matrix(
        as.formula(paste('~', paste(covariates, collapse = " + "))),
        data = pd
    )[, 2:(1 + length(covariates))]

#   Grab the expression for all genes significant in the DEA having significantly
#   associated variants found in the eQTL analysis; convert to long format
express = assays(rse_gene[dea_paired_genes,])$logcounts |>
    cleaningY(mod = mod, P = 1) |>
    as.data.frame() |>
    rownames_to_column("gene_id") |>
    pivot_longer(cols = -gene_id, names_to = "sample_id", values_to = "resid_logcount")

#   Read in genotype data and join with expression data
a = read.plink(paste0(plink_path_prefix, '.bed'))
exp_df = a$genotypes[, dea_paired_variants] |>
    as.data.frame() |>
    rownames_to_column("sample_id") |>
    pivot_longer(
        cols = -sample_id, names_to = "snp_id", values_to = "genotype"
    ) |>
    mutate(
        genotype = factor(3 - as.integer(genotype)),
        gene_id = eqtl$phenotype_id[match(snp_id, eqtl$variant_id)]
    ) |>
    #   Add expression data
    left_join(express, by = c('sample_id', 'gene_id')) |>
    #   Add colData
    left_join(
        colData(rse_gene) |> as_tibble(), by = join_by(sample_id == BrNum)
    ) |>
    filter(sample_id %in% express$sample_id)

#   Plot expression by genotype of each variant with one gene per page and
#   potentially several variants per gene
plot_list_geno = list()
plot_list_dx = list()
plot_list_fraction = list()
for (this_gene in unique(exp_df$gene_id)) {
    this_symbol = rowData(rse_gene)$Symbol[
        match(this_gene, rownames(rse_gene))
    ]

    plot_list_geno[[this_gene]] = exp_df |>
        filter(gene_id == this_gene) |>
        ggplot(mapping = aes(x = genotype, y = resid_logcount)) +
            geom_boxplot(outlier.shape = NA) +
            geom_jitter() +
            facet_wrap(~ snp_id) +
            labs(
                x = "Genotype", y = "Residualized Expression",
                title = this_symbol
            ) +
            theme_bw()
    
    plot_list_dx[[this_gene]] = exp_df |>
        filter(gene_id == this_gene) |>
        #   Don't take duplicate rows if there are multiple SNPs per gene
        filter(snp_id == snp_id[1]) |>
        ggplot(mapping = aes(x = PrimaryDx, y = resid_logcount)) +
            geom_boxplot(outlier.shape = NA) +
            geom_jitter() +
            labs(y = "Residualized Expression", title = this_symbol) +
            theme_bw(base_size = 20)
    
    temp = list()
    for (x_var_name in c("tot.Hb", "tot.Thal")) {
        temp[[x_var_name]] = exp_df |>
            filter(gene_id == this_gene) |>
            ggplot(
                mapping = aes(
                    x = get({{ x_var_name }}), y = resid_logcount,
                    color = genotype
                )
            ) +
            geom_point() +
            geom_smooth(method = lm) +
            theme_bw() +
            labs(
                x = x_var_name, y = "Residualized Expression",
                color = "Genotype"
            )
    }
    plot_list_fraction[[this_gene]] = plot_grid(plotlist = temp, ncol = 2)
}

pdf(file.path(plot_dir, 'expr_by_geno.pdf'))
print(plot_list_geno)
dev.off()

pdf(file.path(plot_dir, 'expr_by_dx.pdf'))
print(plot_list_dx)
dev.off()

pdf(file.path(plot_dir, 'expr_by_geno_fraction.pdf'), width = 10, height = 5)
print(plot_list_fraction)
dev.off()

session_info()
