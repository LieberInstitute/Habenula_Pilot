library("here")
library("sessioninfo")
library("dplyr")
library("purrr")
library("spatialLIBD")
library("ComplexHeatmap")
library("circlize")

out_plot <- here("plots", "12_GSEA")



############## Load objects for gene_set_enrichment_plot_complex ##############

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

###############################################################################


######################### Reproducibility information #########################

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################

