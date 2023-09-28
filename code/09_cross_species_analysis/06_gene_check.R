
library("here")
library("sessioninfo")
library("tidyverse")

load(here("processed-data", 
          "09_cross_species_analysis",
          "Hashikawa_homolog_modeling_results.Rdata"), verbose = TRUE)


search <- list(`Mhb.1` = c("TAC3", "CCK"),
            `Mhb.3` =  c("BHLHE22", "EBF3"),
            `Mhb.2` = c("CHAT"))

map2(search, names(search), function(genes, ct){
  
  hsap_modeling_results$all$enrichment |>
    filter(gene %in% genes) |>
    select(ends_with(ct), ensembl, gene)
  
})

jax_search <- hsap_modeling_results$all$enrichment |>
  filter(gene %in% unlist(search)) |>
  select(ensembl, hs_gene =  gene)

mouse_modeling_results$neuron$enrichment |>
  right_join(jax_search) |>
  select(starts_with("fdr"), ensembl, gene, hs_gene)|>
  pivot_longer(!c(ensembl, gene, hs_gene), names_to = "cell_type", values_to = "fdr", names_prefix = "fdr_") |>
  group_by(gene) |> left_join(mouse_modeling_results$neuron$enrichment |>
  select(starts_with("logFC"), ensembl, gene) |>
  pivot_longer(!c(ensembl, gene), names_to = "cell_type", values_to = "logFC", names_prefix = "logFC_")) |>
  filter(logFC > 0) |>
  arrange(fdr)  |>
  slice(1)

# ensembl gene  hs_gene cell_type    fdr logFC
# <int> <chr> <chr>   <chr>      <dbl> <dbl>
# 1 44789930 Cck   CCK     MHb5      0.113   5.37
# 2 44790260 Chat  CHAT    MHb6      0.0585  4.37
# 3 44791337 Ebf3  EBF3    MHb6      0.263   2.06
# 4 44791974 Tac2  TAC3    MHb6      0.112   3.55
