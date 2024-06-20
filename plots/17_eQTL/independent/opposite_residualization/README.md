In plots in the parent directory, expression is resodualized with two different
models: one for expression vs. diagnosis plots (the model used for DE), and one
for expression vs. SNP genotype plots (the model used for finding eQTLs). This
resulted in unexpected apparently contradictory results when simultaneously
viewing plots for the same information originally residualized differently (e.g. [this](https://github.com/LieberInstitute/Habenula_Pilot/blob/master/plots/17_eQTL/independent/expr_by_dx_eqtls_3_dea_genes_manuscript_FDR1.pdf) and [this](https://github.com/LieberInstitute/Habenula_Pilot/blob/master/plots/17_eQTL/independent/expr_by_geno_eqtls_3_dea_genes_manuscript_FDR1.pdf)
plot). This directory contains plots residualized using the opposite models; i.e.
`expr_by_dx*` plots are residualized by the eQTL model, not DE model, and `expr_by_geno*` plots
are residualized by the DEG model, not eQTL model. These plots were generated interactively
using slightly modified code from [this script](https://github.com/LieberInstitute/Habenula_Pilot/blob/master/code/17_eQTL/04_compare_DEA_GWAS.R).
