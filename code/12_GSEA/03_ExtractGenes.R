library("here")
library("dplyr")
library("purrr")
library("sessioninfo")

out_data <- here("processed-data", "12_GSEA")



###################### Load result tables and sce objects #####################

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
#   enrichTab_FDR05
#   enrichTab_FDR1
#   gene_list

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
# Loading objects:
#   sce_modeling_final_Annotations


## Load gene_set_enrichment_1vsAll_result_tables.rda
load(
    here(
        "processed-data",
        "12_GSEA",
        "gene_set_enrichment_1vsAll_broad_result_tables.rda"
    ),
    verbose = TRUE
)
# Loading objects:
#   enrichTab_br_FDR05
#   enrichTab_br_FDR1
#   gene_list

## Load sce_modeling_broad_Annotations.Rdata
load(
    here(
        "processed-data",
        "05_explore_sce",
        "04_sce_1vALL_modeling",
        "sce_modeling_broad_Annotations.Rdata"
    ),
    verbose = TRUE
)
# Loading objects:
#   sce_modeling_broad_Annotations

###############################################################################



############# Function to extract significant intersection genes ##############

## Extract genes
extract_sig_genes <- function(gsea_sep_tissue, list_set, modeling_results, gene_list, fdr_cut) {
    dims <- dim(modeling_results)[2]
    lapply(gsea_sep_tissue, function(gsea_sep_tissue) {
        gsea_sep_tissue <- cbind(gsea_sep_tissue, modeling_results[, (dims - 1):dims])
        genes <- gsea_sep_tissue %>%
            filter(.[, 3] < fdr_cut & .[, 1] > 0) %>%
            select(ensembl, gene)
        if (list_set == "positive") {
            inter_genes <- intersect(as.vector(genes$ensembl), gene_list$positive)
        } else {
            inter_genes <- intersect(as.vector(genes$ensembl), gene_list$negative)
        }
        genes <- genes %>% filter(ensembl %in% inter_genes)
        return(genes)
    })
}

###############################################################################



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################



