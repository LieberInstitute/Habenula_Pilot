library("here")
library("data.table")
library("dplyr")
library("spatialLIBD")
library("sessioninfo")



##################### Load objects for gene_set_enrichment ####################

## Load sce_modeling_results.Rdata
load(
    here(
        "processed-data",
        "05_explore_sce",
        "04_sce_1vALL_modeling",
        "sce_modeling_results.Rdata"
    ),
    verbose = TRUE
)
# Loading objects:
#   sce_modeling_results

lobstr::obj_size(sce_modeling_results)
# 101.34 MB
class(sce_modeling_results)
# [1] "list"
names(sce_modeling_results)
# [1] "anova"      "enrichment" "pairwise"
dim(sce_modeling_results$enrichment)
# [1] 19799    70

## Load gene list
de_genes <- fread (
    here(
        "processed-data",
        "10_DEA",
        "04_DEA",
        "DEA_Sig-gene_FDR1_qc-totAGene-qSVs-Hb-Thal.tsv"
    ),
    data.table = FALSE,
    stringsAsFactors = FALSE
)

gene_list <- list(
    positive = (de_genes %>%
            filter(logFC > 0) %>%
            select(ensemblID) %>%
            as.vector())$ensemblID,
    negative = (de_genes %>%
            filter(logFC < 0) %>%
            select(ensemblID) %>%
            as.vector())$ensemblID
)

length(gene_list)
# [1] 2
names(gene_list)
# [1] "positive" "negative"
lapply(gene_list, length)
# $positive
# [1] 105
# $negative
# [1] 68

###############################################################################



########################### Run gene_set_enrichment ###########################

enrichTab_FDR05 <- gene_set_enrichment(gene_list = gene_list, modeling_results = sce_modeling_results, model_type = "enrichment", fdr_cut = 0.05)

enrichTab_FDR1 <- gene_set_enrichment(gene_list = gene_list, modeling_results = sce_modeling_results, model_type = "enrichment", fdr_cut = 0.1)

###############################################################################



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################

