library("here")
library("SummarizedExperiment")
library("PCAtools")
library("stringr")
library("ggplot2")
library("cowplot")
library("GGally")
library("ComplexHeatmap")
library("circlize")
library("sessioninfo")

output_path <- here("plots", "10_DEA", "02_DataExploration")



############################# Load rse gene object ############################

load(
    here(
        "processed-data",
        "rse_objects",
        "rse_gene_Habenula_Pilot.rda"
    ),
)

lobstr::obj_size(rse_gene)
# 29.92 MB

class(rse_gene)
# [1] "RangedSummarizedExperiment"
# attr(,"package")
# [1] "SummarizedExperiment"

dim(rse_gene)
# [1] 22756    68

###############################################################################



################### Boxpots: sample_variables - qc_metrics ####################

## NOTE: This section is done with unfiltered and not normalized counts

colData(rse_gene)$library_size <- apply(assay(rse_gene), 2, sum)
colData(rse_gene)$log10_library_size <- log10(colData(rse_gene)$library_size)
colData(rse_gene)$detected_num_genes <- apply(assay(rse_gene), 2, function(x) {
    length(x[which(x > 0)])
})
colData(rse_gene)$abs_ERCCsumLogErr <- abs(colData(rse_gene)$ERCCsumLogErr)

qc_metrics <- c("mitoRate", "rRNA_rate", "overallMapRate", "totalAssignedGene", "concordMapRate", "log10_library_size", "detected_num_genes", "RIN", "abs_ERCCsumLogErr")
sample_variables <- c("PrimaryDx", "Flowcell")


## Function to create boxplots of QC metrics for groups of samples
QC_boxplots <- function(qc_metric, sample_var) {
    if (sample_var == "PrimaryDx") {
        colors <- c("HVYTYBBXX" = "darkmagenta", "HW252BBXX" = "yellow3")
        x_label <- "PrimaryDx"
        sample_var_v2 <- "Flowcell"
    } else if (sample_var == "Flowcell") {
        colors <- c("Schizo" = "darkgoldenrod3", "Control" = "turquoise3")
        x_label <- "Flowcell"
        sample_var_v2 <- "PrimaryDx"
    }

    y_label <- str_replace_all(qc_metric, c("_" = " "))

    data <- data.frame(colData(rse_gene))
    plot <- ggplot(data = data, mapping = aes(x = !!rlang::sym(sample_var), y = !!rlang::sym(qc_metric), color = !!rlang::sym(sample_var_v2))) +
        theme_bw() +
        geom_violin(alpha = 0, size = 0.4, color = "black", width = 0.7) +
        geom_jitter(width = 0.1, alpha = 0.7, size = 2) +
        geom_boxplot(alpha = 0, size = 0.4, width = 0.1, color = "black") +
        scale_color_manual(values = colors) +
        labs(y = y_label, x = x_label) +
        theme(
            axis.title = element_text(size = (9)),
            axis.text = element_text(size = (8))
        )

    return(plot)
}

for (sample_var in sample_variables) {
    i <- 1
    plots <- list()
    for (qc_metric in qc_metrics) {
        plots[[i]] <- QC_boxplots(qc_metric, sample_var)
        i <- i + 1
    }
    plot_grid(plotlist = plots, nrow = 3)

    ggsave(
        paste(output_path, "/QC_boxplots_", sample_var, ".pdf", sep = ""),
        width = 35,
        height = 30,
        units = "cm"
    )
}

###############################################################################



################## Correlation plots: AgeDeath - qc_metrics ###################

## NOTE: This section is done with unfiltered and not normalized counts

rse_gene_df <- data.frame(colData(rse_gene))

## Calculate correlations with Pearson
corrs_age <- sapply(
    qc_metrics,
    function(x) {
        round(cor(rse_gene_df$AgeDeath, rse_gene_df[, x], method = c("pearson")), 3)
    }
)

## Function to plot
plot_correlations <- function(qc_metric) {
    ggplot(rse_gene_df, aes_string(x = "AgeDeath", y = qc_metric)) +
        geom_point(aes_string(colour = "PrimaryDx")) +
        scale_color_manual(values = c("Schizo" = "darkgoldenrod3", "Control" = "turquoise3")) +
        stat_smooth(geom = "line", alpha = 0.7, size = 1.1, span = 0.1, method = lm, show.legend = FALSE) +
        theme_bw() +
        xlab("Age (when death)") +
        ylab(str_replace_all(qc_metric, pattern = "_", replacement = " ")) +
        annotate("text",
            x = max(rse_gene_df$AgeDeath) - max(rse_gene_df$AgeDeath) / 10,
            y = max(rse_gene_df[, qc_metric]) - (max(rse_gene_df[, qc_metric]) - min(rse_gene_df[, qc_metric])) / 10,
            label = paste0("r = ", corrs_age[qc_metric])
        )
}

corr_plots <- lapply(qc_metrics, plot_correlations)
ggsave(
    paste(output_path, "/Corr_AgeDeath_vs_QCmetrics.pdf", sep = ""),
    plot_grid(plotlist = corr_plots, nrow = 3),
    width = 40,
    height = 30,
    units = "cm"
)

###############################################################################



################## PCA and heatmap with colData() variables ###################

## NOTE: This section is done with filtered and normalized counts

## PCA
pca_df <- pca(assays(rse_gene)$logcounts, metadata = colData(rse_gene))

## Plot correlation between PCs and variables (with and without qSVs)
pdf(paste(output_path, "/", "Corr_PCA-Vars-qSVs.pdf", sep = ""), width = 10, height = 10)
eigencorplot(pca_df, metavars = c("PrimaryDx", "AgeDeath", "Flowcell", "mitoRate", "rRNA_rate", "totalAssignedGene", "RIN", "abs_ERCCsumLogErr", "snpPC1", "snpPC2", "snpPC3", "snpPC4", "snpPC5", "tot.Hb", "tot.Thal", "qSV1", "qSV2", "qSV3", "qSV4", "qSV5"))
dev.off()

pdf(paste(output_path, "/", "Corr_PCA-Vars-noqSVs.pdf", sep = ""), width = 10, height = 10)
eigencorplot(pca_df, metavars = c("PrimaryDx", "AgeDeath", "Flowcell", "mitoRate", "rRNA_rate", "totalAssignedGene", "RIN", "abs_ERCCsumLogErr", "snpPC1", "snpPC2", "snpPC3", "snpPC4", "snpPC5", "tot.Hb", "tot.Thal"))
dev.off()

###############################################################################



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
