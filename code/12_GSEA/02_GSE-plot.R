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



#### Function to plot with functions in gene_set_enrichment_plot_complex.R ####

source(
    here(
        "code",
        "12_GSEA",
        "gene_set_enrichment_plot_complex.R"
    )
)

use_gsepc <- function(modeling_results, model_type, gene_list, enrichTab, plot_name) {
    gene_enrichment_count <- get_gene_enrichment_count(model_results = modeling_results, model_type = model_type, bayes_anno = NULL)
    gene_list_count <- get_gene_list_count(gene_list)

    gse_plot <- gene_set_enrichment_plot_complex(
        enrichment = enrichTab,
        gene_count_col = gene_list_count,
        gene_count_row = gene_enrichment_count,
        anno_title_col = "DE Genes",
        anno_title_row = "Cluster\nGenes"
    )

    pdf(paste0(out_plot, "/", plot_name), height = 4, width = 6)
    print(gse_plot)
    dev.off()
}

###############################################################################



################################# Plot results ################################

use_gsepc(
    modeling_results = sce_modeling_final_Annotations,
    model_type = "enrichment",
    gene_list = gene_list,
    enrichTab = enrichTab_FDR05,
    plot_name = "GSEA-1vsAll_FDR05.pdf"
)

use_gsepc(
    modeling_results = sce_modeling_final_Annotations,
    model_type = "enrichment",
    gene_list = gene_list,
    enrichTab = enrichTab_FDR1,
    plot_name = "GSEA-1vsAll_FDR1.pdf"
)

use_gsepc(
    modeling_results = sce_modeling_broad_Annotations,
    model_type = "enrichment",
    gene_list = gene_list,
    enrichTab = enrichTab_br_FDR05,
    plot_name = "GSEA-1vsAll_br_FDR05.pdf"
)

use_gsepc(
    modeling_results = sce_modeling_broad_Annotations,
    model_type = "enrichment",
    gene_list = gene_list,
    enrichTab = enrichTab_br_FDR1,
    plot_name = "GSEA-1vsAll_br_FDR1.pdf"
)

###############################################################################



######################### Reproducibility information #########################

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################

