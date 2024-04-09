library(here)
library(tidyverse)
library(readxl)
library(SummarizedExperiment)
library(snpStats)
library(bigsnpr)
library(sessioninfo)
library(cowplot)
library(data.table)
library(jaffelab)
library(getopt)

#   Read in which tensorQTL run mode is being used
spec <- matrix(
    c("mode", "m", 1, "character", "tensorQTL run mode"),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)

accepted_modes = c('nominal', 'cis', 'independent')
if (!(opt$mode %in% accepted_modes)) {
    stop(
        sprintf(
            "'opt$mode' must be in '%s'.",
            paste(accepted_modes, collapse = "', '")
        )
    )
}

eqtl_path = here(
    'processed-data', '17_eQTL', 'tensorQTL_output', opt$mode, 'FDR05.csv'
)
deg_path = here(
    'processed-data', '10_DEA', '04_DEA',
    'DEA_All-gene_qc-totAGene-qSVs-Hb-Thal.tsv'
)
gwas_narrow_path = here(
    'processed-data', '17_eQTL', 'trubetskoy_gwas_supp_tab1.xls'
)
gwas_narrow_filt_path = here(
    'processed-data', '17_eQTL', 'gwas_narrow_filtered.csv'
)
gwas_wide_path = here(
    "processed-data", "13_MAGMA","GWAS", "scz2022",
    "PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz"
)
gwas_wide_filt_path = here(
    'processed-data', '17_eQTL', 'gwas_wide_filtered.csv'
)
rse_path = here(
    'processed-data', 'rse_objects', 'rse_gene_Habenula_Pilot.rda'
)
plink_path_prefix = here(
    "processed-data", '08_bulk_snpPC', "habenula_genotypes"
)
paired_variants_path = here(
    "processed-data", "17_eQTL", "DEA_paired_variants.txt"
)
raw_geno_path = here(
    'processed-data', '08_bulk_snpPC', 'habenula_genotypes_filt.traw'
)
plot_dir = here('plots', '17_eQTL', opt$mode)

sig_cutoff_deg_explore = c(0.1, 0.05)

#   We can use a stricter cutoff for nominal, where we have more results
if (opt$mode == 'nominal') {
    sig_cutoff_deg_plot = 0.05
} else {
    sig_cutoff_deg_plot = 0.1
}

sig_cutoff_gwas = 5e-8

#   Covariates used in the model for tensorQTL and for DEA, respectively
eqtl_covariates = c(
    "PrimaryDx", 'snpPC1', 'snpPC2', 'snpPC3', 'snpPC4', 'snpPC5'
)
deg_covariates = c(
    'PrimaryDx', 'AgeDeath', 'Flowcell', 'mitoRate', 'rRNA_rate', 'RIN',
    'totalAssignedGene', 'abs_ERCCsumLogErr', 'qSV1', 'qSV2', 'qSV3', 'qSV4',
    'qSV5', 'qSV6', 'qSV7', 'qSV8', 'tot.Hb', 'tot.Thal'
)

lift_over_path = system('which liftOver', intern = TRUE)
dir.create(plot_dir, showWarnings = FALSE)

################################################################################
#   Functions
################################################################################

plot_triad = function(
        rse_gene, dea_paired_variants, mod_deg, mod_eqtl, eqtl, plink, plot_dir,
        plot_prefix, mismatched_snps
    ) {
    #   Grab vector of unique genes paired (via a significant eQTL) with the
    #   variants of interest
    dea_paired_genes = unique(
        eqtl$phenotype_id[
            match(dea_paired_variants, eqtl$variant_id)
        ]
    )

    #   Grab the expression for all genes paired with significant eQTLs;
    #   convert to long format. Residualize using both the tensorQTL and DEA
    #   models
    express_eqtl = assays(rse_gene[dea_paired_genes,])$logcounts |>
        cleaningY(mod = mod_eqtl, P = 1) |>
        as.data.frame() |>
        rownames_to_column("gene_id") |>
        pivot_longer(
            cols = -gene_id, names_to = "sample_id",
            values_to = "resid_logcount_eqtl"
        )
    express_deg = assays(rse_gene[dea_paired_genes,])$logcounts |>
        cleaningY(mod = mod_deg, P = 2) |>
        as.data.frame() |>
        rownames_to_column("gene_id") |>
        pivot_longer(
            cols = -gene_id, names_to = "sample_id",
            values_to = "resid_logcount_deg"
        )

    #   Join genotyping data, expression, and colData
    exp_df = plink$genotypes[, dea_paired_variants] |>
        as.data.frame() |>
        rownames_to_column("sample_id") |>
        pivot_longer(
            cols = -sample_id, names_to = "snp_id", values_to = "genotype"
        ) |>
        mutate(
            genotype = factor(
                #   Ensure 0 is reference, 1 is heterozygous, and 2 is
                #   homozygous for the ALT
                ifelse(
                    snp_id %in% mismatched_snps,
                    as.integer(genotype) - 1,
                    3 - as.integer(genotype)
                )
            ),
            gene_id = eqtl$phenotype_id[match(snp_id, eqtl$variant_id)]
        ) |>
        #   Add expression data residualized from both models
        left_join(express_eqtl, by = c('sample_id', 'gene_id')) |>
        left_join(express_deg, by = c('sample_id', 'gene_id')) |>
        #   Add colData
        left_join(
            colData(rse_gene) |> as_tibble(), by = join_by(sample_id == BrNum)
        ) |>
        filter(sample_id %in% express_eqtl$sample_id)

    #   Plot expression by genotype of each variant with one gene per page and
    #   potentially several variants per gene
    plot_list_geno = list()
    plot_list_dx = list()
    plot_list_fraction = list()
    for (this_gene in unique(exp_df$gene_id)) {
        this_symbol = rowData(rse_gene)$Symbol[
            match(this_gene, rownames(rse_gene))
        ]

        label_df = eqtl |>
            filter(
                phenotype_id == this_gene,
                variant_id %in% exp_df$snp_id
            ) |>
            mutate(sig_label = sprintf(" p = %s", signif(pval_nominal, 3))) |>
            dplyr::rename(snp_id = variant_id)

        plot_list_geno[[this_gene]] = exp_df |>
            filter(gene_id == this_gene) |>
            ggplot() +
                geom_boxplot(
                    mapping = aes(
                        x = genotype, y = resid_logcount_eqtl, color = genotype
                    ),
                    outlier.shape = NA) +
                geom_jitter(
                    mapping = aes(
                        x = genotype, y = resid_logcount_eqtl, color = genotype
                    )
                ) +
                geom_text(
                    data = label_df,
                    mapping = aes(label = sig_label, x = -Inf, y = Inf),
                    hjust = 0,
                    vjust = 1
                ) +
                facet_wrap(~ snp_id) +
                labs(
                    x = "Genotype", y = "Residualized Expression",
                    title = this_symbol
                ) +
                theme_bw() +
                theme(legend.position = "none")
        
        plot_list_dx[[this_gene]] = exp_df |>
            filter(gene_id == this_gene) |>
            #   Don't take duplicate rows if there are multiple SNPs per gene
            filter(snp_id == snp_id[1]) |>
            ggplot(
                    mapping = aes(
                        x = PrimaryDx, y = resid_logcount_deg, color = PrimaryDx
                    )
                ) +
                geom_boxplot(outlier.shape = NA) +
                geom_jitter() +
                labs(y = "Residualized Expression", title = this_symbol) +
                theme_bw(base_size = 20) +
                theme(legend.position = "none")
        
        #   Grab all SNPs associated with this gene
        these_snp_ids = exp_df |>
            filter(gene_id == this_gene) |>
            pull(snp_id) |>
            unique()

        #   Each page will consist of one SNP-gene pair and two plots: one for
        #   habenula and one for thalamus fraction
        for (this_snp_id in these_snp_ids) {
            temp = list()
            this_title = sprintf('%s: %s', this_symbol, this_snp_id)
            for (x_var_name in c("tot.Hb", "tot.Thal")) {
                temp[[x_var_name]] = exp_df |>
                    filter(gene_id == this_gene, snp_id == this_snp_id) |>
                    ggplot(
                        mapping = aes(
                            x = get({{ x_var_name }}), y = resid_logcount_eqtl,
                            color = genotype
                        )
                    ) +
                    geom_point() +
                    geom_smooth(method = lm) +
                    coord_cartesian(xlim = c(0, 1)) +
                    theme_bw(base_size = 20) +
                    labs(
                        x = x_var_name, y = "Residualized Expression",
                        color = "Genotype"
                    )
                
                #   Only want one title and legend, not 2
                if (x_var_name == "tot.Hb") {
                    temp[[x_var_name]] = temp[[x_var_name]] +
                        labs(title = this_title) +
                        theme(legend.position = "none")
                } else {
                    temp[[x_var_name]] = temp[[x_var_name]] +
                        labs(title = " ")
                }
            }
            plot_list_fraction[[this_title]] = plot_grid(plotlist = temp, ncol = 2)
        }
    }

    pdf(file.path(plot_dir, sprintf('expr_by_geno_%s.pdf', plot_prefix)))
    print(plot_list_geno)
    dev.off()

    pdf(file.path(plot_dir, sprintf('expr_by_dx_%s.pdf', plot_prefix)))
    print(plot_list_dx)
    dev.off()

    pdf(
        file.path(
            plot_dir, sprintf('expr_by_geno_fraction_%s.pdf', plot_prefix)
        ),
        width = 14, height = 7
    )
    print(plot_list_fraction)
    dev.off()
}

################################################################################
#   Read in eQTL, DEA, and GWAS results, and the RSE
################################################################################

eqtl = read_csv(eqtl_path, show_col_types = FALSE)

deg_full = read_tsv(deg_path, show_col_types = FALSE)

rse_gene = get(load(rse_path))
colnames(rse_gene) = rse_gene$BrNum
rse_gene$PrimaryDx[rse_gene$PrimaryDx == "Schizo"] = "SCZD"

#   Due to repeated issues with the FTP servers performing the liftover,
#   save filtered copies and load from that rather than lifting over
#   in each interactive session

# gwas_narrow = read_excel(gwas_narrow_path) |>
#     filter(P < sig_cutoff_gwas) |>
#     dplyr::rename(chr = CHR, pos = BP) |>
#     #   Map from hg19 to hg38 and drop anything that fails
#     snp_modifyBuild(lift_over_path, from = 'hg19', to = 'hg38') |>
#     filter(!is.na(pos)) |>
#     #   Construct variant_id from SNP info
#     mutate(
#         variant_id = sprintf(
#             'chr%s:%s:%s:%s',
#             chr,
#             pos,
#             str_split_i(A1A2, '/', 1),
#             str_split_i(A1A2, '/', 2)
#         )
#     )
#
# write_csv(gwas_narrow, gwas_narrow_filt_path)
gwas_narrow = read_csv(gwas_narrow_filt_path, show_col_types = FALSE)

# gwas_wide = fread(gwas_wide_path) |>
#     as_tibble() |>
#     filter(PVAL < sig_cutoff_gwas) |>
#     dplyr::rename(chr = CHROM, pos = POS) |>
#     #   Map from hg19 to hg38 and drop anything that fails
#     snp_modifyBuild(lift_over_path, from = 'hg19', to = 'hg38') |>
#     filter(!is.na(pos)) |>
#     #   Construct variant_id from SNP info
#     mutate(variant_id = sprintf('chr%s:%s:%s:%s', chr, pos, A1, A2))

# write_csv(gwas_wide, gwas_wide_filt_path)
gwas_wide = read_csv(gwas_wide_filt_path, show_col_types = FALSE)

################################################################################
#   Compare each type of result
################################################################################

#-------------------------------------------------------------------------------
#   Info about eQTLS alone
#-------------------------------------------------------------------------------

eqtl_gene = eqtl |>
    group_by(phenotype_id) |>
    summarize(n = n())

#   Significance hardcoded here because it was pre-filtered in 03_adj_p_vals.R
message(
    sprintf("%s significant SNP-gene pairs at FDR<0.05", nrow(eqtl))
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

gwas_paired_genes = list()
for (breadth in c('narrow', 'wide')) {
    if (breadth == 'narrow') {
        gwas = gwas_narrow
    } else {
        gwas = gwas_wide
    }

    message(
        sprintf(
            "Number of %s Trubetskoy et al. GWAS SNPs matching sig eQTLs:",
            breadth
        )
    )
    print(table(gwas$variant_id %in% eqtl$variant_id))

    gwas_paired_genes[[breadth]] = eqtl |>
        filter(variant_id %in% gwas$variant_id) |>
        pull(phenotype_id) |>
        unique()

    message(
        sprintf("Genes paired with %s-GWAS-significant eQTL SNPs:", breadth)
    )
    rowData(rse_gene)$Symbol[
        match(gwas_paired_genes[[breadth]], rownames(rse_gene))
    ] |>
        paste(collapse = ', ') |>
        message()
}


#-------------------------------------------------------------------------------
#   Compare eQTLS with DE genes
#-------------------------------------------------------------------------------

for (sig_cutoff in sig_cutoff_deg_explore) {
    deg =  deg_full |> filter(adj.P.Val < sig_cutoff)
    message(
        sprintf(
            "Number and name of DE genes at FDR<%s implicated in any eQTLs:",
            sig_cutoff
        )
    )
    print(table(deg$gencodeID %in% eqtl_gene$phenotype_id))
    deg$Symbol[deg$gencodeID %in% eqtl_gene$phenotype_id] |>
        paste(collapse = ', ') |>
        message()

    dea_paired_genes = deg$gencodeID[deg$gencodeID %in% eqtl_gene$phenotype_id]

    for (breadth in c('narrow', 'wide')) {
        message(
            sprintf(
                "Genes implicated in the %s GWAS, DEA at FDR<%s, and eQTLs:",
                breadth,
                sig_cutoff
            )
        )

        dea_paired_genes[dea_paired_genes %in% gwas_paired_genes[[breadth]]] |>
            paste(collapse = ', ') |>
            message()
    }
}

################################################################################
#   Read in and preprocess genotype info
################################################################################

#   Read in genotypes
plink = read.plink(paste0(plink_path_prefix, '.bed'))

#   Read in genotype metadata and line up with plink genotypes
geno_raw = read_delim(raw_geno_path, delim = '\t')
geno_raw = geno_raw[match(plink$map$snp.name, geno_raw$SNP),]

#   Keep track of which genotypes should be flipped later, based on code from
#   https://github.com/LieberInstitute/brainseq_phase2/blob/d6779afe6f1509165b4c4eecdfdb7ff7a5d16a19/misc/pull_genotype_data.R#L76-L81
mismatched_snps = plink$map$snp.name[plink$map$allele.2 == geno_raw$COUNTED]

################################################################################
#   Plots exploring how genotype affects expression at select eQTLs
################################################################################

#   Model covariates to later call jaffelab::cleaningY() to residualize
#   expression for downstream plots
pd = as.data.frame(colData(rse_gene))
mod_eqtl <- model.matrix(
        as.formula(paste('~', paste(eqtl_covariates, collapse = " + "))),
        data = pd
    )[, 2:(1 + length(eqtl_covariates))]
mod_deg <- model.matrix(
        as.formula(paste('~', paste(deg_covariates, collapse = " + "))),
        data = pd
    )[, 2:(1 + length(deg_covariates))]

#   After exploring multiple significance cutoffs, choose one for plotting
deg =  deg_full |> filter(adj.P.Val < sig_cutoff_deg_plot)
dea_paired_variants = eqtl |>
    filter(phenotype_id %in% deg$gencodeID) |>
    pull(variant_id)

#   Write paired variants to a text file. This will be read used to subset the
#   big VCF to ensure the below method for reading in genotypes (reading in the
#   plink bed file) works as VariantAnnotation::readVcf() does
if (opt$mode == "nominal") {
    writeLines(dea_paired_variants, paired_variants_path)
}

#   For all modes, plot any SNPs in an eQTL pair paired with a DEG
plot_triad(
    rse_gene = rse_gene,
    dea_paired_variants = dea_paired_variants,
    mod_deg = mod_deg,
    mod_eqtl = mod_eqtl,
    eqtl = eqtl,
    plink = plink,
    plot_dir = plot_dir,
    plot_prefix = sprintf(
        "eqtls_paired_with_dea_genes_FDR%s",
        as.character(sig_cutoff_deg_plot) |> str_split_i('\\.', 2)
    ),
    mismatched_snps
)


if (opt$mode == "independent") {
    #   For independent, plot (11) SNPs overlapping wider GWAS and their paired
    #   genes
    gwas_variants = gwas_wide |>
        filter(variant_id %in% eqtl$variant_id) |>
        pull(variant_id)
    
    plot_triad(
        rse_gene = rse_gene,
        dea_paired_variants = gwas_variants,
        mod_deg = mod_deg,
        mod_eqtl = mod_eqtl,
        eqtl = eqtl,
        plink = plink,
        plot_dir = plot_dir,
        plot_prefix = "wide_gwas_eqtls",
        mismatched_snps
    )

    #   Also plot the top 10 (by significance) eQTLs for independent
    top_eqtls = eqtl |>
        arrange(FDR) |>
        slice_head(n = 10) |>
        pull(variant_id)
    
    plot_triad(
        rse_gene = rse_gene,
        dea_paired_variants = top_eqtls,
        mod_deg = mod_deg,
        mod_eqtl = mod_eqtl,
        eqtl = eqtl,
        plink = plink,
        plot_dir = plot_dir,
        plot_prefix = "top_10_eqtls",
        mismatched_snps
    )
}

session_info()
