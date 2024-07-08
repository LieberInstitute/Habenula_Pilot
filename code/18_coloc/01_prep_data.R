library(here)
library(tidyverse)
library(SummarizedExperiment)
library(sessioninfo)
library(data.table)
library(bigsnpr)
library(miniparquet)


bim_path = here(
   'processed-data', '08_bulk_snpPC', 'v3', 'habenula_R.9_MAF.05.RSann.bim'
)
rse_path = here(
    'processed-data', 'rse_objects', 'rse_gene_Habenula_Pilot.rda'
)
gwas_wide_path = here(
    "processed-data", "13_MAGMA","GWAS", "scz2022",
    "PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz"
)
eqtl_dir = here("processed-data", "17_eQTL", "tensorQTL_output", "nominal")
window_radius = 500000

lift_over_path = system('which liftOver', intern = TRUE)

gwas_wide = fread(gwas_wide_path) |>
    as_tibble() |>
    dplyr::rename(chr = CHROM, pos = POS) |>
    #   Map from hg19 to hg38 and drop anything that fails
    snp_modifyBuild(lift_over_path, from = 'hg19', to = 'hg38') |>
    filter(!is.na(pos)) |>
    #   Construct variant_id from SNP info
    mutate(variant_id = sprintf('chr%s:%s:%s:%s', chr, pos, A1, A2)) |>
    #   Only care about SNPs also evaluated in eQTL data
    filter(variant_id %in% eqtl$variant_id)

parquet_files <- list.files(
    eqtl_dir,
    pattern = "\\.parquet$",
    full.names = TRUE
)

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
