library("SummarizedExperiment")
library("sessioninfo")
library("tidyverse")
library("VariantAnnotation")
library("jaffelab")
library("here")
library("recount")
library("data.table")

source(here("eqtl", "code", "rse_to_bed.R"))

rse_path = here(
    'processed-data', 'rse_objects', 'rse_gene_Habenula_Pilot.rda'
)
pca_path = here(
    'processed-data', '03_bulk_pca', '02_multiregion_PCA',
    'Multi_region_PCs.Rdata'
)
expected_covariates = c(
    "PrimaryDx", 'snpPC1', 'snpPC2', 'snpPC3', 'snpPC4', 'snpPC5'
)
out_dir = here("processed-data", "17_eQTL", "tensorQTL_input")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

rse = get(load(rse_path, verbose = TRUE))

samples_split <- map(rse_split, colnames)

#### Covariate Data ####
load(pca_path, verbose = TRUE)
pcs = colData(rse) |>
    as_tibble() |>
    dplyr::select(RNum) |>
    left_join(pca_tab, by = "RNum") |>
    dplyr::select(RNum, starts_with("PC")) |>
    column_to_rownames("RNum")

corner(pcs)

message(Sys.time(), " - Format covariates")
## Phenodata
pd = as.data.frame(colData(rse)[, expected_covariates])

pd <- model.matrix(
        as.formula(paste('~', paste(expected_covariates, collapse = " + "))),
        data = pd
    )[, 2:(1 + length(expected_covariates))]

covars = cbind(pd, pcs) |>
    t() |>
    as.data.frame() |>
    rownames_to_column("id")
corner(covars)

write_tsv(covars, file = file.path(out_dir, "covariates.txt"))


#### Expression Data ####

## test TSS fix 8/24/23
# bed <- rse_to_bed(rse_split$amyg)
# corner(bed)
# tail(colnames(bed))
# bed |> filter(ID == "ENSG00000000003.14") |> dplyr::select(ID, `#Chr`, `start`, `end`)
# ID #Chr     start       end
# 1 ENSG00000000003.14 chrX 100639991 100639992

message(Sys.time(), " - logcounts to bed")
expression_fn <- map(features, function(feat) map(regions, ~ here("eqtl", "data", "tensorQTL_input", "expression_bed", paste0(feat, "_", .x, ".bed.gz"))))

expression_bed <- map2(list(rse, rse_exon, rse_jxn, rse_tx), features, function(rse, feat) {
    message(Sys.time(), " - ", feat)
    rse_split <- map(regions, ~ rse[, rse$BrainRegion == .x])
    expr_bed <- map(rse_split, rse_to_bed)
    return(expr_bed)
})

## double check the output
map_depth(expression_bed, 2, corner)

message(Sys.time(), " - Write bed tables to .bed.gz files")
walk2(expression_bed, expression_fn, function(expr, fn) {
    walk2(expr, fn, ~ data.table::fwrite(.x, .y,
        sep = "\t",
        quote = FALSE, row.names = FALSE
    ))
})


#### VCF ####
# message(Sys.time(), " - Split VCF")
# risk_vcf <- readVcf(here("eqtl", "data", "risk_snps", "LIBD_maf01_gwas_BPD.vcf.gz"))
# risk_vcf
# risk_vcf_split <- map(rse_split, ~ risk_vcf[, .x$genoSample])
# map(risk_vcf_split, dim)
# 
# vcf_fn <- map(regions, ~ here("eqtl", "data", "risk_snps", paste0("LIBD_maf01_gwas_BPD_", .x, ".vcf.gz")))
# walk2(risk_vcf_split, vcf_fn, ~ writeVcf(.x, .y))
# 
# ## plink commands
# map(vcf_fn, ~ paste("plink --make-bed --output-chr chrM --vcf", .x, "--out", gsub(".vcf.gz", "", .x)))
# 

## check

# map2(bed, risk_vcf_split, ~ all(colnames(.x[, 5:ncol(.x)]) == colnames(.y)))
# map2(bed, covars, ~ all(colnames(.x[, 5:ncol(.x)]) == colnames(.y[[1]][, 2:ncol(.y[[1]])])))


#### prep interaction csv ####
message(Sys.time(), " - Prep interaction data")
walk2(rse_split, regions, function(rse, region) {
    cell_fractions <- colData(rse)[, c("Astro", "Endo", "Macro", "Micro", "Mural", "Oligo", "OPC", "Tcell", "Excit", "Inhib")]
    cell_fractions <- as.data.frame(cell_fractions)
    rownames(cell_fractions) <- rse$genoSample
    write.csv(cell_fractions, file = here("eqtl", "data", "tensorQTL_input", "interaction", paste0("cell_fraction_", region, ".csv")))
})

## create shell commands ##
# sgejobs::job_single("tensorqtl_risk_snps",
#     create_shell = TRUE, queue = "bluejay", memory = "50G",
#     command = "python tensorqtl_risk_snps.py"
# )

# sgejobs::job_single('01_convert_rdata_to_tensorqtl_format', create_shell = TRUE, memory = '50G', command = "Rscript 01_convert_rdata_to_tensorqtl_format.R")

## Reproducibility information
print("Reproducibility information:")
options(width = 120)
session_info()