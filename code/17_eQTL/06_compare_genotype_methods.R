#   This script was run interactively to compare methods for reading in
#   genotypes: from plink bed vs. VariantAnnotation::readVcf()

library(here)
library(tidyverse)
library(SummarizedExperiment)
library(snpStats)
library(VariantAnnotation)

plink_path = here(
    "processed-data", '08_bulk_snpPC', "habenula_genotypes.bed"
)
paired_variants_path = here(
    'processed-data', '17_eQTL', 'DEA_paired_variants.txt'
)
small_vcf_path = here(
    'processed-data', '17_eQTL', 'paired_DEA_eQTL_SNPs.vcf'
)
#   We'll also read just the MAF in for all SNPs to know for the manuscript
#   how many SNPs remain after the MAF>0.05 cutoff
big_vcf_path = here('processed-data', '08_bulk_snpPC', 'habenula_genotypes.vcf.gz')

paired_variants = readLines(paired_variants_path)

#   Read in genotypes from the plink bed and tidy
a = read.plink(plink_path)
a = a$genotypes[, paired_variants] |>
    as.data.frame() |>
    rownames_to_column("sample_id") |>
    pivot_longer(
        cols = -sample_id, names_to = "snp_id", values_to = "genotype"
    ) |>
    mutate(genotype = as.integer(genotype)) |>
    as_tibble()

#   Read in genotypes from the subsetted VCF and tidy
b = geno(readVcf(small_vcf_path))$GT |>
    as.data.frame() |>
    rownames_to_column("snp_id") |>
    pivot_longer(
        cols = -snp_id, names_to = "sample_id", values_to = "genotype"
    ) |>
    mutate(
        genotype = case_when(
            genotype == "0|0" ~ 0,
            genotype == "1|1" ~ 2,
            TRUE ~ 1
        )
    ) |>
    as_tibble()

#   Merge tibbles so genotype values are separate columns that can be compared
#   side by side
combined_df = left_join(
    a |> dplyr::rename(genotype_a = genotype),
    b |> dplyr::rename(genotype_b = genotype),
    by = c('sample_id', 'snp_id')
)

#   None are equal as-is
table(combined_df$genotype_a == combined_df$genotype_b)

#   Are the numbers just in reverse order?
combined_df_temp = combined_df |>
    mutate(genotype_a = 3 - genotype_a)

#   Ok, that was easy
all(combined_df_temp$genotype_a == combined_df_temp$genotype_b)

#   Read in MAF for all SNPs to have the number passing the cutoff for the
#   manuscript
a = readInfo(big_vcf_path, x = 'MAF')
message(
    sprintf(
        "%s of %s SNPs remain after MAF>0.05 cutoff applied.",
        sum(a > 0.05),
        length(a)
    )
)
# 6763560 of 12393872 SNPs remain after MAF>0.05 cutoff applied.
