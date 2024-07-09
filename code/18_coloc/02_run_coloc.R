library(here)
library(tidyverse)
library(sessioninfo)
library(coloc)

in_dir = here("processed-data", "18_coloc", "temp_inputs")
out_dir = here("processed-data", "18_coloc", "coloc_individual_out")
gene_path = here("processed-data", "18_coloc", "all_genes.txt")

dir.create(out_dir, showWarnings = FALSE)

gene = system(
    paste('awk "NR == $SLURM_ARRAY_TASK_ID"', gene_path), intern = TRUE
)
message("Proceeding with gene ", gene)

dataset_list = readRDS(file.path(in_dir, paste0(gene, '.rds')))

results = coloc.abf(dataset_list[['eqtl']], dataset_list[['gwas']])

#   Grab the 95% credible set of causal SNPs
o <- order(results$results$SNP.PP.H4, decreasing=TRUE)
cs <- cumsum(results$results$SNP.PP.H4[o])
w <- which(cs > 0.95)[1]
cred_set = results$results[o,][1:w, c('snp', 'SNP.PP.H4')]

#   Write SNPs and the probabilities both traits are associated and share the
#   SNPs as causal (if any such SNPs exist)
if (nrow(cred_set > 0)) {
    write_csv(cred_set, file.path(out_dir, paste0(gene, '.csv')))
} else {
    message("No SNPs in the 95% credible set.")
}
