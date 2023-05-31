library("here")
library("SummarizedExperiment")
library("edgeR")
library("limma")
library("sessioninfo")



############################# Load rse gene object ############################

load(
    here(
        "preprocessed_data",
        "count_data_bukola",
        "rse_gene_Roche_Habenula_qcAndAnnotated_n69.Rdata"
    ),
    verbose = TRUE
)

lobstr::obj_size(rse_gene)
# 40.63 MB

class(rse_gene)
# [1] "RangedSummarizedExperiment"
# attr(,"package")
# [1] "SummarizedExperiment"

rse_gene
# class: RangedSummarizedExperiment
# dim: 58037 69
# metadata(0):
# assays(1): counts
# rownames(58037): ENSG00000223972.5 ENSG00000227232.5 ...
#   ENSG00000210195.2 ENSG00000210196.2
# rowData names(10): Length gencodeID ... NumTx gencodeTx
# colnames: NULL
# colData names(68): RNum RIN ... Flowcell hasGenotype

dim(colData(rse_gene))
# [1] 69 68

colData(rse_gene)$library_size <- apply(assay(rse_gene), 2, sum)
colData(rse_gene)$log10_library_size <- log10(colData(rse_gene)$library_size)
colData(rse_gene)$detected_num_genes <- apply(assay(rse_gene), 2, function(x) {
    length(x[which(x > 0)])
})

###############################################################################



############################## Gene normalization #############################

length(which(assay(rse_gene) == 0)) * 100 / (dim(rse_gene)[1] * dim(rse_gene)[2])
# [1] 45.30061

# hist(assay(rse_gene), breaks = 50)

assays(rse_gene, withDimnames = FALSE)$logcounts <- edgeR::cpm(calcNormFactors(rse_gene, method = "TMM"), log = TRUE, prior.count = 0.5)

# hist(assays(rse_gene)$logcounts, breaks = 50, prob = TRUE)

###############################################################################
######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################


