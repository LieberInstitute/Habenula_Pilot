library("here")
library("SummarizedExperiment")
library("PCAtools")
library("stringr")
library("ggplot2")
library("cowplot")
library("ComplexHeatmap")
library("circlize")
library("sessioninfo")

output_path <- here("plots", "10_DEA", "02_DataExploration")



############################ Load rse gene objects ############################

load(
    here(
        "processed-data",
        "rse_objects",
        "rse_gene_Habenula_Pilot.rda"
    ),
    verbose = TRUE
)
# Loading objects:
#   rse_gene

lobstr::obj_size(rse_gene)
# 29.92 MB

class(rse_gene)
# [1] "RangedSummarizedExperiment"
# attr(,"package")
# [1] "SummarizedExperiment"

rse_gene
# class: RangedSummarizedExperiment
# dim: 22756 68
# metadata(0):
# assays(2): counts logcounts
# rownames(22756): ENSG00000227232.5 ENSG00000278267.1 ...
#   ENSG00000210195.2 ENSG00000210196.2
# rowData names(11): Length gencodeID ... gencodeTx MGI_Symbol
# colnames(68): R18346 R18347 ... R18422 R18423
# colData names(109): BrNum RNum ... snpPC9 snpPC10

dim(colData(rse_gene))
# [1] 68 109

names(colData(rse_gene))
#   [1] "BrNum"                          "RNum"
#   [3] "RIN"                            "Brain.Region"
#   [5] "AgeDeath"                       "Sex"
#   [7] "Race"                           "PrimaryDx"
#   [9] "FQCbasicStats"                  "perBaseQual"
#  [11] "perTileQual"                    "perSeqQual"
#  [13] "perBaseContent"                 "GCcontent"
#  [15] "Ncontent"                       "SeqLengthDist"
#  [17] "SeqDuplication"                 "OverrepSeqs"
#  [19] "AdapterContent"                 "KmerContent"
#  [21] "SeqLength_R1"                   "percentGC_R1"
#  [23] "phred20.21_R1"                  "phred48.49_R1"
#  [25] "phred76.77_R1"                  "phred100.101_R1"
#  [27] "phredGT30_R1"                   "phredGT35_R1"
#  [29] "Adapter50.51_R1"                "Adapter70.71_R1"
#  [31] "Adapter88.89_R1"                "SeqLength_R2"
#  [33] "percentGC_R2"                   "phred20.21_R2"
#  [35] "phred48.49_R2"                  "phred76.77_R2"
#  [37] "phred100.101_R2"                "phredGT30_R2"
#  [39] "phredGT35_R2"                   "Adapter50.51_R2"
#  [41] "Adapter70.71_R2"                "Adapter88.89_R2"
#  [43] "ERCCsumLogErr"                  "bamFile"
#  [45] "trimmed"                        "numReads"
#  [47] "numMapped"                      "numUnmapped"
#  [49] "overallMapRate"                 "concordMapRate"
#  [51] "totalMapped"                    "mitoMapped"
#  [53] "mitoRate"                       "totalAssignedGene"
#  [55] "gene_Assigned"                  "gene_Unassigned_Ambiguity"
#  [57] "gene_Unassigned_MultiMapping"   "gene_Unassigned_NoFeatures"
#  [59] "gene_Unassigned_Unmapped"       "gene_Unassigned_MappingQuality"
#  [61] "gene_Unassigned_FragmentLength" "gene_Unassigned_Chimera"
#  [63] "gene_Unassigned_Secondary"      "gene_Unassigned_Nonjunction"
#  [65] "gene_Unassigned_Duplicate"      "rRNA_rate"
#  [67] "Flowcell"                       "hasGenotype"
#  [69] "sum"                            "detected"
#  [71] "subsets_Mito_sum"               "subsets_Mito_detected"
#  [73] "subsets_Mito_percent"           "subsets_Ribo_sum"
#  [75] "subsets_Ribo_detected"          "subsets_Ribo_percent"
#  [77] "library_size"                   "log10_library_size"
#  [79] "detected_num_genes"             "abs_ERCCsumLogErr"
#  [81] "Astrocyte"                      "Endo"
#  [83] "Excit.Thal"                     "Inhib.Thal"
#  [85] "LHb"                            "MHb"
#  [87] "Microglia"                      "Oligo"
#  [89] "OPC"                            "tot.Hb"
#  [91] "tot.Thal"                       "qSV1"
#  [93] "qSV2"                           "qSV3"
#  [95] "qSV4"                           "qSV5"
#  [97] "qSV6"                           "qSV7"
#  [99] "qSV8"                           "snpPC1"
# [101] "snpPC2"                         "snpPC3"
# [103] "snpPC4"                         "snpPC5"
# [105] "snpPC6"                         "snpPC7"
# [107] "snpPC8"                         "snpPC9"
# [109] "snpPC10"

unique(colData(rse_gene)$PrimaryDx)
# [1] "Schizo"  "Control"

table(colData(rse_gene)$PrimaryDx)
# Control  Schizo
#      33      35

table(colData(rse_gene)$Sex)
#  M
# 68

table(colData(rse_gene)$Race)
# CAUC
#   68

table(colData(rse_gene)$Flowcell)
# HVYTYBBXX HW252BBXX
#        34        34

###############################################################################



################### Boxpots: sample_variables - qc_metrics ####################

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
eigencorplot(
    pca_df,
    metavars = c("PrimaryDx", "AgeDeath", "Flowcell", "mitoRate", "rRNA_rate", "totalAssignedGene", "RIN", "abs_ERCCsumLogErr", "snpPC1", "snpPC2", "snpPC3", "snpPC4", "snpPC5", "tot.Hb", "tot.Thal", "qSV1", "qSV2", "qSV3", "qSV4", "qSV5"),
    col = c("#005f73", "#0a9396", "#94d2bd", "#e9d8a6", "#ee9b00", "#ca6702", "#bb3e03")
)
dev.off()

pdf(paste(output_path, "/", "Corr_PCA-Vars-noqSVs.pdf", sep = ""), width = 10, height = 10)
eigencorplot(
    pca_df,
    metavars = c("PrimaryDx", "AgeDeath", "Flowcell", "mitoRate", "rRNA_rate", "totalAssignedGene", "RIN", "abs_ERCCsumLogErr", "snpPC1", "snpPC2", "snpPC3", "snpPC4", "snpPC5", "tot.Hb", "tot.Thal"),
    col = c("#005f73", "#0a9396", "#94d2bd", "#e9d8a6", "#ee9b00", "#ca6702", "#bb3e03")
)
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
