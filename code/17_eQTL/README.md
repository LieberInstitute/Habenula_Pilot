This README describes results from the `04_compare_DEA_GWAS.*` scripts, which explore `tensorQTL` eQTL results, and check for overlap with [our DEGs](https://github.com/LieberInstitute/Habenula_Pilot/blob/d3810a75cb9dba46b14e6cb51a3bdac44c8abe1d/code/17_eQTL/04_compare_DEA_GWAS.R#L33-L36) and PGC3 GWAS SNPs at [two levels](https://github.com/LieberInstitute/Habenula_Pilot/blob/d3810a75cb9dba46b14e6cb51a3bdac44c8abe1d/code/17_eQTL/04_compare_DEA_GWAS.R#L37-L43): a larger set of ~20k SNPs, and a smaller set of just 313.

## Overview of results

This table summarizes the number of significant results when comparing overlap of eQTLs with DEGs and GWAS SNPs. While [3 plots are produced for each eQTL set](#plot-organization-and-naming-conventions), links to the first plot (against diagnosis) are provided in the table where they exist.

| eQTL type | eQTL FDR | num sig genes total | num sig genes overlap FDR < 0.05 DEGs | num sig genes overlap FDR < 0.1 DEGs | num sig SNP-gene pairs | num sig SNPs overlap GWAS 313 p < 5e-8 | num sig SNPs overlap GWAS 20k p < 5e-8 | trifectas at DEG FDR < 0.05 | trifectas at DEG FDR < 0.1 |
| ----- | ----- | ------- | ------- | ------- | ------- | ------- | ------- | ------- | ------- |
| nominal | 0.05 | 3163 | [5](https://github.com/LieberInstitute/Habenula_Pilot/blob/master/plots/17_eQTL/nominal/expr_by_dx_eqtls_paired_with_dea_genes_FDR05.pdf) | 23 | 147744 | 7 | 4285 | 0 | 0 |
| nominal | 0.01 | 1311 | 1 | 12 | 86656 | 3 | 3354 | 0 | 0 |
| cis | 0.05 | 728 | 0 | 3 | 728 | 0 | 11 | 0 | 0 |
| independent | 0.05 | 728 | 0 | [3](https://github.com/LieberInstitute/Habenula_Pilot/blob/master/plots/17_eQTL/independent/expr_by_dx_eqtls_paired_with_dea_genes_FDR1.pdf) | 742 | 0 | [11](https://github.com/LieberInstitute/Habenula_Pilot/blob/master/plots/17_eQTL/independent/expr_by_dx_wide_gwas_eqtls.pdf) | 0 | 0 |

## Plot organization and naming conventions

### Prefix

For a given set of eQTLs, there are three plots under `plots/17_eQTL/{eQTL level}`:

1. Starting with `expr_by_dx_`: for each gene paired with any eQTL in the set, plot residualized expression against diagnosis.
2. Starting with `expr_by_geno_`: for each eQTL in the set, sort by paired gene and plot residualized expression against genotype
3. Starting with `expr_by_geno_fraction_`: for each eQTL in the set, plot residualized expression against habenula and thalamus fraction, colored by genotype

### Base name

The set of eQTLs to plot determine the remainder of the filename; e.g. `expr_by_dx_eqtls_paired_with_dea_genes_FDR05.pdf` is one of three plots for any eQTL paired with a gene found to be differentially expressed at FDR < 0.05.

## Statistical significance of results

There are generally 3 significance thresholds used in plots (one for eQTLs, one for DEGs, and one for GWAS SNPs):

- **eQTLs**: always FDR < 0.05
- **DEGs**: FDR cutoff given in the filename if applicable (e.g. `expr_by_dx_eqtls_paired_with_dea_genes_FDR1.pdf` uses FDR < 0.1)
- **GWAS SNPs**: always the genome-wide p-value cutoff of 5e-8 used in the PGC3 manuscript