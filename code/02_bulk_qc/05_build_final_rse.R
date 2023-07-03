library("here")
library("SummarizedExperiment")
library("sessioninfo")



############################# Load rse gene object ############################

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

## Add some metrics to colData
colData(rse_gene)$library_size <- apply(assay(rse_gene), 2, sum)
colData(rse_gene)$log10_library_size <- log10(colData(rse_gene)$library_size)
colData(rse_gene)$detected_num_genes <- apply(assay(rse_gene), 2, function(x) {
    length(x[which(x > 0)])
})
colData(rse_gene)$abs_ERCCsumLogErr <- abs(colData(rse_gene)$ERCCsumLogErr)

## Add total proportion of Hb (LHb + MHb) and total Thal (Excit.Thal + Inhib.Thal)
est_prop <- t(est_prop$bulk.props)

est_prop <- cbind(est_prop, est_prop[, colnames(est_prop) == "LHb"] + est_prop[, colnames(est_prop) == "MHb"])
est_prop <- cbind(est_prop, est_prop[, colnames(est_prop) == "Excit.Thal"] + est_prop[, colnames(est_prop) == "Inhib.Thal"])

colnames(est_prop)[10:11] <- c("tot.Hb", "tot.Thal")

## Delete sample Br5572/R18424
rse_gene <- rse_gene[, rse_gene$BrNum != "Br5572"]
est_prop <- est_prop[rownames(est_prop) != "R18424", ]

###############################################################################



########### Add deconvolution, SNP PCs and qSVa results to colData ############

## Add cell proportions and qSVAs
rownames(est_prop) == rownames(colData(rse_gene)) ## I'm checking if samples are in the same order
rownames(qSVAs) == rownames(colData(rse_gene))

colData(rse_gene) <- cbind(colData(rse_gene), est_prop, qSVAs)

## Add SNP PCs
colData(rse_gene) <- merge(colData(rse_gene), as.data.frame(snpPCs), by = "BrNum")
colnames(rse_gene) <- colData(rse_gene)$RNum ## The merge deletes the colnames so I'm adding them back

###############################################################################



######################## Save rse object with qsva data #######################

save(rse_gene,
    file = here(
        "processed-data",
        "rse_objects",
        "rse_gene_Habenula_Pilot.rda"
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

