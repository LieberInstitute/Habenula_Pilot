library("here")
library("SummarizedExperiment")
library("sessioninfo")

out_data <- here("processed-data", "rse_objects")



#################### Load rse gene, tx, jx and exon objects ###################

load(
    here(
        "processed-data",
        "02_bulk_qc",
        "count_data_bukola",
        "rse_gene_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata"
    ),
    verbose = TRUE
)
# Loading objects:
#   rse_gene

load(
    here(
        "processed-data",
        "02_bulk_qc",
        "count_data_bukola",
        "rse_tx_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata"
    ),
    verbose = TRUE
)
# Loading objects:
#   rse_tx

load(
    here(
        "processed-data",
        "02_bulk_qc",
        "count_data_bukola",
        "rse_jx_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata"
    ),
    verbose = TRUE
)
# Loading objects:
#   rse_jx

load(
    here(
        "processed-data",
        "02_bulk_qc",
        "count_data_bukola",
        "rse_exon_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata"
    ),
    verbose = TRUE
)
# Loading objects:
#   rse_exon

###############################################################################



################# Load deconvolution, SNP PCs and qSVa results ################

## Load deconvolution results
load(
    here(
        "processed-data",
        "06_deconvolution",
        "run_Bisque",
        "est_prop_split_Hb_annotations.RDATA"
    ),
    verbose = TRUE
)
# Loading objects:
#   est_prop


## Load SNP PCs
snpPCs <- read.table(
    here(
        "processed-data",
        "08_bulk_snpPC",
        "habenula_genotypes_filt.snpPCs.tab"
    ),
    header = TRUE
)

## Load qSVa
qSVAs <- read.table(
    here(
        "processed-data",
        "11_bulk_qSVA",
        "qSVA.tsv"
    ),
    header = TRUE
)

###############################################################################



####################### Prepare data to add to the model ######################

## Add some metrics only to the rse_gene colData
colData(rse_gene)$library_size <- apply(assay(rse_gene), 2, sum)
colData(rse_gene)$log10_library_size <- log10(colData(rse_gene)$library_size)
colData(rse_gene)$detected_num_genes <- apply(assay(rse_gene), 2, function(x) {
    length(x[which(x > 0)])
})

## Add the absolute value of ERCCsumLogErr to all the rses
colData(rse_gene)$abs_ERCCsumLogErr <- abs(colData(rse_gene)$ERCCsumLogErr)
colData(rse_tx)$abs_ERCCsumLogErr <- abs(colData(rse_gene)$ERCCsumLogErr)
colData(rse_exon)$abs_ERCCsumLogErr <- abs(colData(rse_gene)$ERCCsumLogErr)
colData(rse_jx)$abs_ERCCsumLogErr <- abs(colData(rse_gene)$ERCCsumLogErr)

## Add total proportion of Hb (LHb + MHb) and total Thal (Excit.Thal + Inhib.Thal)
est_prop <- t(est_prop$bulk.props)
est_prop <- cbind(est_prop, est_prop[, colnames(est_prop) == "LHb"] + est_prop[, colnames(est_prop) == "MHb"])
est_prop <- cbind(est_prop, est_prop[, colnames(est_prop) == "Excit.Thal"] + est_prop[, colnames(est_prop) == "Inhib.Thal"])

colnames(est_prop)[10:11] <- c("tot.Hb", "tot.Thal")

## Delete sample Br5572/R18424 from all rse objects and deconvolution data
rse_gene <- rse_gene[, rse_gene$BrNum != "Br5572"]
rse_tx <- rse_tx[, rse_tx$BrNum != "Br5572"]
rse_exon <- rse_exon[, rse_exon$BrNum != "Br5572"]
rse_jx <- rse_jx[, rse_jx$BrNum != "Br5572"]

est_prop <- est_prop[rownames(est_prop) != "R18424", ]

###############################################################################



########### Add deconvolution, SNP PCs and qSVa results to colData ############

add_vars <- function(rse) {
    ## Add cell proportions and qSVAs
    stopifnot(identical(rownames(est_prop), rownames(colData(rse)))) ## Check if samples are in the same order
    stopifnot(identical(rownames(qSVAs), rownames(colData(rse))))

    colData(rse) <- cbind(colData(rse), est_prop, qSVAs)

    ## Add SNP PCs
    merged_col <- merge(colData(rse), as.data.frame(snpPCs), by = "BrNum", sort = FALSE)
    stopifnot(identical(rse$RNum, merged_col$RNum))
    colData(rse) <- merged_col
    colnames(rse) <- colData(rse)$RNum ## The merge deletes the colnames so I'm adding them back

    return(rse)
}

rse_gene <- add_vars(rse_gene)
colnames(rse_tx) <- colData(rse_tx)$RNum
rse_tx <- add_vars(rse_tx)
rse_exon <- add_vars(rse_exon)
rse_jx <- add_vars(rse_jx)

###############################################################################



################# Save rse objects with qsva and SNP PCs data #################

save_rse <- function(rse, rse_name) {
    rse_orName <- deparse(substitute(rse))
    print(rse_orName)
    assign(rse_orName, rse)
    save(
        list = rse_orName,
        file = paste0(out_data, "/", rse_orName, "_Habenula_Pilot.rda")
    )
}

save_rse(rse_gene, "rse_gene")
save_rse(rse_tx, "rse_tx")
save_rse(rse_exon, "rse_exon")
save_rse(rse_jx, "rse_jx")

###############################################################################



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
