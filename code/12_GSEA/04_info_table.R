library(tidyverse)
library(data.table)
library(here)
library(sessioninfo)

## Load gene list
de_genes <- fread(
    here(
        "processed-data",
        "10_DEA",
        "04_DEA",
        "DEA_Sig-gene_FDR1_qc-totAGene-qSVs-Hb-Thal.tsv"
    ),
    data.table = FALSE,
    stringsAsFactors = FALSE
)

head(de_genes)

de_genes2 <- de_genes |>
    filter(adj.P.Val < 0.1) |>
    select(ensemblID, DGE_logFC = logFC, DGE_t = t, DGE_fdr = adj.P.Val) |>
    as_tibble() |>
    mutate(DGE_direction = ifelse(DGE_t > 0, "positive", "negative"))

 ## matches
de_genes2 |> count(DGE_direction)
# direction     n
# <chr>     <int>
# 1 negative     68
# 2 positive    105


## Load gene_set_enrichment_1vsAll_result_tables.rda
load(
    here(
        "processed-data",
        "12_GSEA",
        "gene_set_enrichment_1vsAll_result_tables.rda"
    ),
    verbose = TRUE
)

# Loading objects:
# enrichTab_FDR05
# enrichTab_FDR1
# gene_list

enrichTab_FDR1 |> filter(ID == "all", NumSig != 0)
#    OR       Pval       test NumSig SetSize  ID model_type fdr_cut
# 1  1.1885593 0.32907199  Astrocyte     12     128 all enrichment     0.1
# 2  1.4834351 0.15668629       Endo     10     128 all enrichment     0.1
# 3  1.8173905 0.08644021 Excit.Thal      8     128 all enrichment     0.1
# 4  1.7871165 0.02801306 Inhib.Thal     16     128 all enrichment     0.1
# 5 19.3378892 0.05671375      LHb.3      1     128 all enrichment     0.1
# 6  2.1736418 0.37364082      MHb.2      1     128 all enrichment     0.1
# 7  1.3518216 0.14876470  Microglia     18     128 all enrichment     0.1
# 8  0.3759059 0.96763630      Oligo      2     128 all enrichment     0.1
# 9  1.5061022 0.22137544        OPC      6     128 all enrichment     0.1

## Load sce_modeling_final_Annotations.Rdata
load(
    here(
        "processed-data",
        "05_explore_sce",
        "04_sce_1vALL_modeling",
        "sce_modeling_final_Annotations.Rdata"
    ),
    verbose = TRUE
)

ct_enrich <- sce_modeling_final_Annotations$enrichment |>
    select(starts_with("fdr"), ensemblID = ensembl, gene) |>
    pivot_longer(!c(gene, ensemblID), names_to ="cell_type", values_to = "ct_enrichment_fdr", names_prefix = "fdr_") |>
    left_join(sce_modeling_final_Annotations$enrichment |>
                  select(starts_with("t_stat_"), ensemblID = ensembl) |>
                  pivot_longer(!ensemblID, names_to ="cell_type", values_to = "ct_enrichment_t", names_prefix = "t_stat_"),
              relationship = "many-to-many") |>
    filter(ct_enrichment_fdr < 0.1, ct_enrichment_t > 0)

ct_enrich |> count(cell_type)
# cell_type      n
# <chr>      <int>
# 1 Astrocyte   1587
# 2 Endo        1073
# 3 Excit.Thal   704
# 4 Inhib.Thal  1472
# 5 LHb.1         37
# 6 LHb.2         32
# 7 LHb.3          9
# 8 LHb.5          1
# 9 LHb.6         19
# 10 LHb.7          6
# 11 MHb.1         60
# 12 MHb.2         72
# 13 MHb.3         45
# 14 Microglia   2142
# 15 OPC          628
# 16 Oligo        799

gse_info <- de_genes2 |>
    inner_join(ct_enrich) |>
    arrange(cell_type, ct_enrichment_fdr)

gse_info |> count(cell_type)
# cell_type      n
# <chr>      <int>
# 1 Astrocyte     12
# 2 Endo          10
# 3 Excit.Thal     8
# 4 Inhib.Thal    16
# 5 LHb.3          1
# 6 MHb.2          1
# 7 Microglia     18
# 8 OPC            6
# 9 Oligo          2


write.csv(gse_info, file = here("processed-data", "12_GSEA", "GSEA_FDR1_all_info.csv"))

