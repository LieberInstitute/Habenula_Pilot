library(here)
library(tidyverse)
library(sessioninfo)

results_path = here("processed-data", "18_coloc", "coloc_results.rds")
out_path = here("processed-data", "18_coloc", "filtered_results.csv")
h4_cutoff = 0.8

results = readRDS(results_path)

shared_causal_variant = sapply(
    results, function(x) x$summary[['PP.H4.abf']] > h4_cutoff
)
message(sprintf("Genes with a shared causal variant (cutoff > %s):", h4_cutoff))
table(shared_causal_variant)

all_genes = names(results)

#   Add gene to each tibble of results
results_df = lapply(
    1:length(results),
    function(i) {
        as_tibble(results[[i]]$results) |>
            mutate(gene = all_genes[i])
    }
)

do.call(rbind, results_df) |>
    #   Only take genes with sufficient probability of shared causal variant
    filter(gene %in% all_genes[shared_causal_variant]) |>
    #   For each gene, take just the set of SNPs in the 95% credible set
    group_by(gene) |>
    arrange(desc(SNP.PP.H4)) |>
    mutate(h4_cumsum = cumsum(SNP.PP.H4)) |>
    filter(h4_cumsum <= 0.95) |>
    ungroup() |>
    #   Save just the gene, SNP, and shared-causal-variant probability
    select(gene, snp, SNP.PP.H4) |>
    write_csv(out_path)

session_info()
