library("tidyverse")
library("jaffelab")
library("miniparquet")
library("sessioninfo")
library("here")
library("qs")
library("getopt")

#   Read in which tensorQTL run mode is being used
spec <- matrix(
    c("mode", "m", 1, "character", "tensorQTL run mode"),
    c("covariate", "c", 2, "character", "covariate applicable for interaction model"),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)

accepted_modes = c('nominal', 'cis', 'independent', 'interaction')
if (!(opt$mode %in% accepted_modes)) {
    stop(
        sprintf(
            "'opt$mode' must be in '%s'.",
            paste(accepted_modes, collapse = "', '")
        )
    )
}

if (opt$mode == "interaction") {
    out_dir_suffix = sprintf('%s_%s', opt$mode, opt$covariate)
} else {
    out_dir_suffix = opt$mode
}

out_dir = here("processed-data", "17_eQTL", "tensorQTL_output", out_dir_suffix)
p_val_cutoff = 0.05

if (opt$mode == 'nominal') {
    parquet_files <- list.files(
        out_dir,
        pattern = "\\.parquet$",
        full.names = TRUE
    )

    eqtl_out <- do.call("rbind", map(parquet_files, parquet_read)) |>
        mutate(FDR = p.adjust(pval_nominal, "fdr"))
} else if (opt$mode == "interaction") {
    parquet_files <- list.files(
        out_dir,
        pattern = "\\.parquet$",
        full.names = TRUE
    )

    eqtl_out <- do.call("rbind", map(parquet_files, parquet_read)) |>
        mutate(FDR = p.adjust(pval_gi, "fdr"))
} else if (opt$mode == 'cis') {
    eqtl_out = read_csv(
            file.path(out_dir, paste0(opt$mode, '_out.csv')),
            show_col_types = FALSE
        ) |>
        mutate(FDR =  p.adjust(pval_beta, "fdr"))
} else if (opt$mode == "independent") {
    eqtl_out = read_csv(
            file.path(out_dir, paste0(opt$mode, '_out.csv')),
            show_col_types = FALSE
        ) |>
        mutate(FDR =  p.adjust(pval_perm, "fdr"))
} else { 
    stop(sprintf("Unsupported 'opt$mode' = '%s'.", opt$mode))
}

message("n pairs: ", nrow(eqtl_out))

# filter
eqtl_out_filtered <- eqtl_out |>
    filter(FDR < p_val_cutoff)
message("n pairs FDR<", p_val_cutoff, ": ", nrow(eqtl_out_filtered))

fn <- here(
    out_dir,
    paste0(
        "FDR", as.character(p_val_cutoff) |> str_extract('\\.(.*)$', group = 1)
    )
)

#   Save as CSV and the faster qs
write_csv(eqtl_out_filtered, file = paste0(fn, ".csv"))
qsave(eqtl_out_filtered, file = paste0(fn, ".qs"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
