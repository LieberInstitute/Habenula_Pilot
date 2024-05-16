library("SummarizedExperiment")
library("sessioninfo")
library("tidyverse")
library("jaffelab")
library("here")
library("data.table")

rse_path = here(
    'processed-data', 'rse_objects', 'rse_gene_Habenula_Pilot.rda'
)
pca_path = here('processed-data', '03_bulk_pca', 'PCs.rds')
snp_pcs_path = here(
    'processed-data', '08_bulk_snpPC',
    'Hb_gt_merged_R.9_MAF.05_ann_filt.snpPCs.tab'
)
expected_covariates = c(
    "PrimaryDx", 'snpPC1', 'snpPC2', 'snpPC3', 'snpPC4', 'snpPC5',
    paste0('PC', 1:12)
)
out_dir = here("processed-data", "17_eQTL", "tensorQTL_input")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
#   Function definitions
################################################################################

# Convert assay from a SummarizedExperiment to tensorQTL compatible bed file
rse_to_bed <- function(rse, assay_name = "logcounts") {
    stopifnot(assay_name %in% assayNames(rse))

    rr_df <- as.data.frame(rowRanges(rse)) |>
        tibble::rownames_to_column("ID") |>
        dplyr::mutate(
            start = ifelse(strand == "+", start, end),
            end = start + 1
        ) |>
        dplyr::arrange(seqnames, start) |>
        dplyr::select(`#Chr` = seqnames, start, end, ID)
  
    message(
        Sys.time(), " - Converting assay '", assay_name, "' to bed format: (",
        nrow(rse), ", ", ncol(rse),")"
    )
    counts <- SummarizedExperiment::assays(rse)[[assay_name]]

    counts <- rr_df |> 
        dplyr::left_join(
            as.data.frame(counts) |> 
                tibble::rownames_to_column("ID"), 
                by = "ID"
            )
    return(counts)
}

################################################################################
#   Compute and write 'covariates.txt'
################################################################################

rse = get(load(rse_path, verbose = TRUE))

colData(rse) = colData(rse) |>
    as_tibble() |>
    #   Overwrite SNP PCs with the most recent values computed from the properly
    #   filtered genotyping data
    select(!matches('^snpPC')) |>
    left_join(
        read_tsv(snp_pcs_path, show_col_types = FALSE) |>
            #   Fix the name of one donor
            mutate(BrNum = ifelse(BrNum == "Br0983", "Br983", BrNum)),
        by = 'BrNum'
    ) |>
    #   Add in gene PCs
    left_join(
        readRDS(pca_path) |> as.data.frame() |> rownames_to_column("RNum"),
        by = "RNum"
    ) |>
    DataFrame()
colnames(rse) = rse$BrNum

#   Also write colData to CSV, to have easy access to potential interaction
#   covariates
col_data = as_tibble(colData(rse))
colnames(col_data) = gsub('\\.', '_', colnames(col_data)) 
write_csv(col_data, file.path(out_dir, 'colData.csv'))

message(Sys.time(), " - Format covariates")

#   Model covariates and save
pd = as.data.frame(colData(rse)[, expected_covariates])
covars <- model.matrix(
        as.formula(paste('~', paste(expected_covariates, collapse = " + "))),
        data = pd
    )[, 2:(1 + length(expected_covariates))] |>
    t() |>
    as.data.frame() |>
    rownames_to_column("id")
corner(covars)

write_tsv(covars, file = file.path(out_dir, "covariates.txt"))

################################################################################
#   Write logcounts to a ".bed.gz" file
################################################################################

message(Sys.time(), " - Format logcounts as bed in memory")
expression_bed <- rse_to_bed(rse)
corner(expression_bed)

message(Sys.time(), " - Write bed table to .bed.gz file")
data.table::fwrite(
    expression_bed, here(out_dir, "logcounts.bed.gz"), sep = "\t",
    quote = FALSE, row.names = FALSE
)

## Reproducibility information
print("Reproducibility information:")
options(width = 120)
session_info()
