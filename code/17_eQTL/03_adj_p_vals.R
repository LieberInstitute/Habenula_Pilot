library("tidyverse")
library("jaffelab")
library("miniparquet")
library("sessioninfo")
library("here")
library("qs")

out_dir = here("processed-data", "17_eQTL", "tensorQTL_output")
p_val_cutoff = 0.05

parquet_files <- list.files(
    out_dir,
    pattern = "", # TODO
    full.names = TRUE
)

eqtl_out <- do.call("rbind", map(parquet_files, parquet_read)) |>
    mutate(FDR = p.adjust(pval_nominal, "fdr")) 
message("n pairs: ", nrow(eqtl_out))

# filter
eqtl_out_filtered <- eqtl_out |>
    filter(FDR < cutoff)
message("n pairs FDR<", cutoff, ": ", nrow(eqtl_out_filtered))
    
fn <- here(out_dir, "FDR05")

#   Save as CSV and the faster qs
write_csv(eqtl_out_filtered, file = paste0(fn, ".csv"))
qsave(eqtl_out_filtered, file = paste0(fn, ".qs"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
