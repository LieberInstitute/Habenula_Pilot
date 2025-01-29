#   After the first round of manuscript reviews, we were asked to perform
#   threshold-free comparisons of DE results between habenula and other
#   brain regions. This script performs these comparisons

library(here)
library(tidyverse)
library(BiocFileCache)
library(readxl)
library(sessioninfo)

hab_de_path = here(
    "processed-data", "10_DEA", "04_DEA",
    'DEA_All-gene_qc-totAGene-qSVs-Hb-Thal.tsv'
)
bsp2_dlpfc_de_path = here(
    "raw-data", "15_compare_vs_BrainSeq",
    "dxStats_dlpfc_filtered_qSVA_geneLevel_noHGoldQSV_matchDLPFC.rda"
)
bsp2_hippo_de_path = here(
    "raw-data", "15_compare_vs_BrainSeq",
    "dxStats_hippo_filtered_qSVA_geneLevel_noHGoldQSV_matchHIPPO.rda"
)
bsp3_caudate_de_url = "https://caudate-eqtl.s3.us-west-2.amazonaws.com/BrainSeq_Phase3_Caudate_DifferentialExpression_DxSZ_all.txt.gz"
dg_de_path = here(
    "raw-data", "15_compare_vs_BrainSeq", "41593_2020_604_MOESM2_ESM.xlsx"
)

################################################################################
#   Load and clean DE results for all brain regions
################################################################################

#   Habenula
hab_de = read_tsv(hab_de_path, show_col_types = FALSE) |>
    select(ensemblID, t, P.Value, adj.P.Val)

#   BSP2 DLPFC
load(bsp2_dlpfc_de_path, verbose = TRUE)
bsp2_dlpfc_de = outGene |>
    as_tibble() |>
    select(ensemblID, t, P.Value, adj.P.Val)
rm(outGene0, outGeneNoAdj, outGene)

#   BSP2 Hippocampus
load(bsp2_hippo_de_path, verbose = TRUE)
bsp2_hippo_de = outGene |>
    as_tibble() |>
    select(ensemblID, t, P.Value, adj.P.Val)
rm(outGene0, outGeneNoAdj, outGene)

#   BSP3 Caudate
bfc <- BiocFileCache()
bsp3_caudate_file <- BiocFileCache::bfcrpath(
    bfc, bsp3_caudate_de_url, exact = TRUE
)
bsp3_caudate_de = read_tsv(bsp3_caudate_file) |>
    filter(Type == "Gene") |>
    select(ensemblID, t, P.Value, adj.P.Val)

#   DG
dg = read_xlsx(dg_de_path, sheet = "Table_S10") |>
    select(ensemblID, SZ_t, SZ_P.Value, SZ_adj.P.Val) |>
    rename(t = SZ_t, P.Value = SZ_P.Value, adj.P.Val = SZ_adj.P.Val)