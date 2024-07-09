library(here)
library(tidyverse)
library(sessioninfo)
library(data.table)
library(bigsnpr)
library(miniparquet)
library(coloc)
library(BiocParallel)

gwas_wide_path = here(
    "processed-data", "13_MAGMA","GWAS", "scz2022",
    "PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz"
)
eqtl_dir = here("processed-data", "17_eQTL", "tensorQTL_output", "nominal")
out_path = here("processed-data", "18_coloc", "coloc_results.rds")

lift_over_path = system('which liftOver', intern = TRUE)
n_cores = as.integer(Sys.getenv('SLURM_CPUS_PER_TASK'))

################################################################################
#   Prepare GWAS and eQTL data
################################################################################

#   Read in full set of GWAS SNPs and lift to hg38
message(Sys.time(), ' | Reading in GWAS SNPs and lifting to hg38...')
gwas_wide = fread(gwas_wide_path) |>
    as_tibble() |>
    dplyr::rename(chr = CHROM, pos = POS) |>
    #   Map from hg19 to hg38 and drop anything that fails
    snp_modifyBuild(lift_over_path, from = 'hg19', to = 'hg38') |>
    filter(!is.na(pos)) |>
    #   Construct variant_id from SNP info
    mutate(variant_id = sprintf('chr%s:%s:%s:%s', chr, pos, A1, A2))

#   Some SNPs have a different number of case and control donors as evidence,
#   while coloc requires a single value across all SNPs. Here we take the
#   median, noting that the variance is extremely low (i.e. most SNPs hover
#   around the median value)
prop_cases_vec = gwas_wide$NCAS / (gwas_wide$NCAS + gwas_wide$NCON)
prop_cases = median(prop_cases_vec)
stopifnot(sd(prop_cases_vec) < 1e-2)

parquet_files <- list.files(
    eqtl_dir,
    pattern = "\\.parquet$",
    full.names = TRUE
)
message(Sys.time(), ' | Reading in eQTL results and merging with GWAS...')
coloc_df <- do.call("rbind", map(parquet_files, parquet_read)) |>
    as_tibble() |>
    #   Only consider eQTLs whose SNPs are in GWAS data
    filter(variant_id %in% gwas_wide$variant_id) |>
    #   Reformat for input to coloc
    dplyr::rename(snp = variant_id, beta_eqtl = slope) |>
    mutate(varbeta_eqtl = slope_se ** 2) |>
    select(phenotype_id, snp, beta_eqtl, varbeta_eqtl) |>
    #   Grab beta and varbeta values from the GWAS data at the same SNPs
    left_join(
        gwas_wide |>
            dplyr::rename(snp = variant_id, beta_gwas = BETA) |>
            mutate(varbeta_gwas = SE ** 2) |>
            select(snp, beta_gwas, varbeta_gwas),
        by = 'snp'
    )

################################################################################
#   Run coloc and save results in individual RDS files per gene
################################################################################

run_coloc = function(gene, coloc_df) {
    this_coloc_df = coloc_df |>
        filter(phenotype_id == gene)

    eqtl = this_coloc_df |>
        dplyr::rename(beta = beta_eqtl, varbeta = varbeta_eqtl) |>
        select(snp, beta, varbeta) |>
        as.list()
    eqtl$type = "quant"
    eqtl$sdY = 1

    gwas = this_coloc_df |>
        dplyr::rename(beta = beta_gwas, varbeta = varbeta_gwas) |>
        select(snp, beta, varbeta) |>
        as.list()
    gwas$type = "cc"
    gwas$s = prop_cases

    return(coloc.abf(eqtl, gwas))
}

message(Sys.time(), ' | Running coloc...')

all_genes = unique(coloc_df$phenotype_id)

results_list = bplapply(
    all_genes,
    run_coloc,
    coloc_df = coloc_df,
    BPPARAM = MulticoreParam(n_cores)
)
names(results_list) = all_genes

saveRDS(results_list, out_path)

session_info()
