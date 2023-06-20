library("here")
library("SummarizedExperiment")
library("ggplot2")
library("cowplot")
library("GGally")
library("sessioninfo")

output_path <- here("plots", "08_bulk_snpPC", "02_explore_snpPCs")

############################ Load rse gene objects ############################

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

lobstr::obj_size(rse_gene)
# 30.17 MB

class(rse_gene)
# [1] "RangedSummarizedExperiment"
# attr(,"package")
# [1] "SummarizedExperiment"

snpPCs <- read.table(
    here(
        "processed-data",
        "08_bulk_snpPC",
        "habenula_genotypes_n69.snpPCs.tab"
    ),
    header = TRUE
)

lobstr::obj_size(snpPCs)
# 11.80 kB

class(snpPCs)
# [1] "data.frame"

###############################################################################



################################ Plotting PCs #################################

colData(rse_gene)$library_size <- apply(assay(rse_gene), 2, sum)
colData(rse_gene)$log10_library_size <- log10(colData(rse_gene)$library_size)
colData(rse_gene)$detected_num_genes <- apply(assay(rse_gene), 2, function(x) {
    length(x[which(x > 0)])
})
colData(rse_gene)$abs_ERCCsumLogErr <- abs(colData(rse_gene)$ERCCsumLogErr)

cols_rse <- colData(rse_gene)[order(colData(rse_gene)$BrNum), ]
cols_rse <- cols_rse[cols_rse$BrNum != "Br6158",]
snpPCs <- snpPCs[order(snpPCs$BrNum),]

cols_rse <- cbind(snpPCs[,-1], cols_rse)
cols_rse <- as.data.frame(cols_rse)

sample_variables <- c("PrimaryDx", "Flowcell", "AgeDeath")

pca_samplevars <- function(sample_var) {
    pair_pc <- ggpairs(cols_rse,
        columns = 1:4, aes_string(color = sample_var, alpha = 0.5),
        upper = list(continuous = wrap("points", size = 2))
    )

    if (sample_var == "Flowcell") {
        for (i in 1:4) {
            for (j in 1:4) {
                if (j != i) {
                    plot <- ggplot(cols_rse, aes_string(x = paste0("snpPC", i), y = paste0("snpPC", j), color = sample_var, label = "BrNum")) +
                        geom_point() +
                        scale_color_manual(values = c("HVYTYBBXX" = "darkmagenta", "HW252BBXX" = "yellow3")) +
                        geom_text() +
                        theme_bw()
                    pair_pc[i, j] <- plot
                }
            }
        }
    } else if (sample_var == "PrimaryDx") {
        for (i in 1:4) {
            for (j in 1:4) {
                if (j != i) {
                    plot <- ggplot(cols_rse, aes_string(x = paste0("snpPC", i), y = paste0("snpPC", j), color = sample_var, label = "BrNum")) +
                        geom_point() +
                        scale_color_manual(values = c("Schizo" = "darkgoldenrod3", "Control" = "turquoise3")) +
                        geom_text() +
                        theme_bw()
                    pair_pc[i, j] <- plot
                }
            }
        }
    } else if (sample_var == "AgeDeath") {
        for (i in 1:4) {
            for (j in 1:4) {
                if (j != i) {
                    plot <- ggplot(cols_rse, aes_string(x = paste0("snpPC", i), y = paste0("snpPC", j), color = sample_var, label = "BrNum")) +
                        geom_point() +
                        scale_color_gradient(low = "#F19E93", high = "#00433F") +
                        geom_text() +
                        theme_bw()
                    pair_pc[i, j] <- plot
                }
            }
        }
    }

    ggsave(
        paste(output_path, "/PCA_", sample_var, ".pdf", sep = ""),
        pair_pc,
        width = 30,
        height = 34,
        units = "cm"
    )
}

lapply(sample_variables, pca_samplevars)

## Matrix of PCs

qc_metrics <- c("mitoRate", "rRNA_rate", "overallMapRate", "totalAssignedGene", "concordMapRate", "log10_library_size", "detected_num_genes", "RIN", "abs_ERCCsumLogErr")

pca_qcmetrics <- function(qc_metric) {
    pair_pc <- ggpairs(cols_rse,
        columns = 1:4, aes_string(color = qc_metric, alpha = 0.5),
        upper = list(continuous = wrap("points", size = 2))
    )

    for (i in 1:4) {
        for (j in 1:4) {
            if (j != i) {
                plot <- ggplot(cols_rse, aes_string(x = paste0("snpPC", i), y = paste0("snpPC", j), color = qc_metric, label = "BrNum")) +
                    geom_point() +
                    scale_color_gradient(low = "yellow", high = "darkblue") +
                    geom_text() +
                    theme_bw()
                pair_pc[i, j] <- plot
            }
        }
    }

    ggsave(
        paste(output_path, "/PCA_", qc_metric, ".pdf", sep = ""),
        pair_pc,
        width = 30,
        height = 34,
        units = "cm"
    )
}

lapply(qc_metrics, pca_qcmetrics)

###############################################################################



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
