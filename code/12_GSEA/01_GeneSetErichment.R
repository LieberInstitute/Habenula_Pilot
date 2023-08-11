library("here")
library("data.table")
library("dplyr")
library("spatialLIBD")
library("sessioninfo")

out_data <- here("processed-data", "12_GSEA")
if (!dir.exists(out_data)) dir.create(out_data)



##################### Load objects for gene_set_enrichment ####################

message(Sys.time(), " - loading objects")

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

lobstr::obj_size(sce_modeling_final_Annotations)
# 101.34 MB
class(sce_modeling_final_Annotations)
# [1] "list"
names(sce_modeling_final_Annotations)
# [1] "anova"      "enrichment" "pairwise"
dim(sce_modeling_final_Annotations$enrichment)
# [1] 19799    70

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

lobstr::obj_size(sce_modeling_broad_Annotations)
# 32.07 MB
class(sce_modeling_broad_Annotations)
# [1] "list"
names(sce_modeling_broad_Annotations)
# [1] "anova"      "enrichment" "pairwise"
dim(sce_modeling_broad_Annotations$enrichment)
# [1] 19324    38


## Load marker_stats_top_25_genes.Rdata
load(
    here(
        "processed-data",
        "06_deconvolution",
        "run_Bisque",
        "marker_stats_top_25_genes.Rdata"
    ),
    verbose = TRUE
)
# Loading objects:
#   marker_stats

lobstr::obj_size(marker_stats)
# 5.23 MB
class(marker_stats)
# [1] "tbl_df"     "tbl"        "data.frame"
dim(marker_stats)
# [1] 23082    18

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

gene_list <- list(
    all = de_genes$ensemblID,
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
# [1] 3
names(gene_list)
# [1] "all"      "positive" "negative"
lapply(gene_list, length)
# $all
# [1] 173
# $positive
# [1] 105
# $negative
# [1] 68

###############################################################################



################## Run gene_set_enrichment (broad and final) ##################

message(Sys.time(), " - running gene_set_enrichment() for broad and final annotations")

enrichTab_FDR05 <- gene_set_enrichment(gene_list = gene_list, modeling_results = sce_modeling_final_Annotations, model_type = "enrichment", fdr_cut = 0.05)

enrichTab_FDR1 <- gene_set_enrichment(gene_list = gene_list, modeling_results = sce_modeling_final_Annotations, model_type = "enrichment", fdr_cut = 0.1)

enrichTab_br_FDR05 <- gene_set_enrichment(gene_list = gene_list, modeling_results = sce_modeling_broad_Annotations, model_type = "enrichment", fdr_cut = 0.05)

enrichTab_br_FDR1 <- gene_set_enrichment(gene_list = gene_list, modeling_results = sce_modeling_broad_Annotations, model_type = "enrichment", fdr_cut = 0.1)

###############################################################################



################# Run gene_set_enrichment for 25 marker genes #################

## Since the object marker_stats doesn't have the same structure as a normal
## modelig_results object, I'm using other function to run it.Which was created
## by Louise Huuki

message(Sys.time(), " - running marker_gene_set_enrichment() for top 25 marker genes")

source(
    here(
        "code",
        "12_GSEA",
        "marker_gene_set_enrichment.R"
    )
)

enrichTab_mrk <- marker_gene_set_enrichment(gene_list = gene_list, marker_stats = marker_stats)

###############################################################################



################### Save gene_set_enrichment results to rda ###################

message(Sys.time(), " - saving all results")

save(enrichTab_FDR05,
    enrichTab_FDR1,
    gene_list,
    file = here(
        out_data,
        "gene_set_enrichment_1vsAll_result_tables.rda"
    )
)

save(enrichTab_br_FDR05,
    enrichTab_br_FDR1,
    gene_list,
    file = here(
        out_data,
        "gene_set_enrichment_1vsAll_broad_result_tables.rda"
    )
)

save(enrichTab_mrk,
    gene_list,
    file = here(
        out_data,
        "gene_set_enrichment_1vsAll_25mrkr_result_tables.rda"
    )
)

###############################################################################



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
