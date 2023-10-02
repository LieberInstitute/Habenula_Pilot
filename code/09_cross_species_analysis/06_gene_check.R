
library("here")
library("sessioninfo")
library("tidyverse")
library("SingleCellExperiment")

load(here("processed-data", 
          "09_cross_species_analysis",
          "Hashikawa_homolog_modeling_results.Rdata"), verbose = TRUE)

# load(file = here("processed-data","sce_objects", "sce_Habenula_Pilot.Rdata"), verbose = TRUE)


search <- list(`Mhb.1` = c("TAC3", "CCK"),
            `Mhb.3` =  c("BHLHE22", "EBF3"),
            `Mhb.2` = c("CHAT"),
            `Mhb` = "CHRNB4",
            `Lhb.3` = "ONECUT2",
            `LHb.3` = "MCOLN3",
            `LHb.1` = "SEMA3D",
            `LHb` = "HTR4",
            `LHb.6` = "ESRP1",
            `LHb.4` = "TLE2",
            `LHb.2` = "CRH",
            `LHb` = "HTR4")

map2(search, names(search), function(genes, ct){
  
  hsap_modeling_results$all$enrichment |>
    filter(gene %in% genes) |>
    select(ends_with(ct), ensembl, gene)
  
})

## ensembl is really JAXid 
search_tb <- tibble(cell_type_search = names(search), gene = search) |>
  unnest(gene) |>
  left_join(hsap_modeling_results$all$enrichment) |>
  filter(gene %in% unlist(search)) |>
  select(cell_type_search, ensembl, hs_gene =  gene)

# cell_type_search  ensembl hs_gene
# <chr>               <int> <chr>  
#   1 Mhb.1            44791974 TAC3   
# 2 Mhb.1            44789930 CCK    
# 3 Mhb.3            44787249 BHLHE22
# 4 Mhb.3            44791337 EBF3   
# 5 Mhb.2            44790260 CHAT   
# 6 Lhb3             44800479 ONECUT2
# 7 LHb3                   NA MCOLN3 
# 8 LHb1             44784746 SEMA3D 
# 9 LHb              44795328 HTR4   
# 10 LHb6             44794749 ESRP1  
# 11 LHb2                   NA CRH    
# 12 LHb              44795328 HTR4  


search_tb |> filter(!ensembl %in% mouse_modeling_results$neuron$enrichment$ensembl)
search_tb |> filter(!ensembl %in% mouse_modeling_results$all$enrichment$ensembl)

enrich_long_mm_neuron <- mouse_modeling_results$neuron$enrichment |>
  select(starts_with("fdr"), ensembl, gene)|>
  pivot_longer(!c(ensembl, gene), names_to = "cell_type", values_to = "fdr", names_prefix = "fdr_") |>
  left_join(mouse_modeling_results$neuron$enrichment |>
              select(starts_with("logFC"), ensembl) |>
              pivot_longer(!c(ensembl), names_to = "cell_type", values_to = "logFC", names_prefix = "logFC_")) 

enrich_long_mm_all <- mouse_modeling_results$all$enrichment |>
  select(starts_with("fdr"), ensembl, gene)|>
  pivot_longer(!c(ensembl, gene), names_to = "cell_type", values_to = "fdr", names_prefix = "fdr_") |>
  left_join(mouse_modeling_results$all$enrichment |>
              select(starts_with("logFC"), ensembl) |>
              pivot_longer(!c(ensembl), names_to = "cell_type", values_to = "logFC", names_prefix = "logFC_")) 

## search of genes of interest

search_tb |> 
  left_join(enrich_long_mm_neuron, relationship = "many-to-many") |>
  filter(logFC > 1) |>
  group_by(gene) |>
  arrange(fdr) |>
  dplyr::slice(1)

# cell_type_search  ensembl hs_gene gene    cell_type    fdr logFC
# <chr>               <int> <chr>   <chr>   <chr>      <dbl> <dbl>
#   1 Mhb.1            44789930 CCK     Cck     MHb5      0.113   5.37
# 2 Mhb.2            44790260 CHAT    Chat    MHb6      0.0585  4.37
# 3 Mhb.3            44791337 EBF3    Ebf3    MHb6      0.263   2.06
# 4 LHb              44795328 HTR4    Htr4    MHb2      0.981   1.98
# 5 Lhb.3            44800479 ONECUT2 Onecut2 LHb2      0.171   2.67
# 6 LHb.1            44784746 SEMA3D  Sema3d  LHb4      0.417   2.12
# 7 Mhb.1            44791974 TAC3    Tac2    MHb6      0.112   3.55
# 8 LHb.4            44789237 TLE2    Tle2    LHb2      0.120   2.46

search_tb |> 
  left_join(enrich_long_mm_all, relationship = "many-to-many") |>
  filter(logFC > 1) |>
  group_by(gene) |>
  arrange(fdr) |>
  dplyr::slice(1)

# cell_type_search  ensembl hs_gene gene    cell_type         fdr logFC
# <chr>               <int> <chr>   <chr>   <chr>           <dbl> <dbl>
#   1 Mhb.3            44787249 BHLHE22 Bhlhe22 OPC3      0.123        2.88
# 2 Mhb.1            44789930 CCK     Cck     Neuron7   0.493        4.16
# 3 Mhb.2            44790260 CHAT    Chat    Neuron7   0.194        4.41
# 4 Mhb.3            44791337 EBF3    Ebf3    Microglia 0.227        2.09
# 5 LHb.6            44794749 ESRP1   Esrp1   OPC1      0.000000233  5.69
# 6 LHb              44795328 HTR4    Htr4    Neuron7   0.496        2.39
# 7 Lhb.3            44800479 ONECUT2 Onecut2 Oligo3    0.171        4.39
# 8 LHb.1            44784746 SEMA3D  Sema3d  OPC3      0.0759       5.13
# 9 Mhb.1            44791974 TAC3    Tac2    Neuron7   0.278        5.26
# 10 LHb.4            44789237 TLE2    Tle2    OPC2      0.489        2.95

enrich_long_hs_all <- hsap_modeling_results$all$enrichment |>
  select(starts_with("fdr"), ensembl, gene)|>
  pivot_longer(!c(ensembl, gene), names_to = "cell_type", values_to = "fdr", names_prefix = "fdr_") |>
  left_join(hsap_modeling_results$all$enrichment |>
              select(starts_with("logFC"), ensembl) |>
              pivot_longer(!c(ensembl), names_to = "cell_type", values_to = "logFC", names_prefix = "logFC_")) 


search_tb |> 
  left_join(enrich_long_hs_all, relationship = "many-to-many") |>
  filter(logFC > 1) |>
  group_by(gene) |>
  arrange(fdr) |>
  dplyr::slice(1)

# cell_type_search  ensembl hs_gene gene    cell_type           fdr logFC
# <chr>               <int> <chr>   <chr>   <chr>             <dbl> <dbl>
# 1 Mhb.3            44787249 BHLHE22 BHLHE22 MHb.3      0.000171      6.59
# 2 Mhb.1            44789930 CCK     CCK     Excit.Thal 0.108         3.60
# 3 Mhb.2            44790260 CHAT    CHAT    MHb.2      0.000214      5.38
# 4 Mhb.3            44791337 EBF3    EBF3    MHb.3      0.171         7.19
# 5 LHb.6            44794749 ESRP1   ESRP1   LHb.6      0.368         4.48
# 6 LHb              44795328 HTR4    HTR4    LHb.7      0.801         2.29
# 7 Lhb.3            44800479 ONECUT2 ONECUT2 LHb.1      0.178         4.64
# 8 LHb.1            44784746 SEMA3D  SEMA3D  Astrocyte  0.203         3.49
# 9 Mhb.1            44791974 TAC3    TAC3    MHb.1      0.0000000200  7.46
# 10 LHb.4            44789237 TLE2    TLE2    Inhib.Thal 0.112         3.10
