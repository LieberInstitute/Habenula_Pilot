library("here")
library("SummarizedExperiment")
library("recount")
library("dplyr")
library("ggplot2")
library("cowplot")
library("stringr")
library("RColorBrewer")
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

names(colData(rse_gene))
#  [1] "RNum"                           "RIN"
#  [3] "Brain.Region"                   "BrNum"
#  [5] "AgeDeath"                       "Sex"
#  [7] "Race"                           "PrimaryDx"
#  [9] "FQCbasicStats"                  "perBaseQual"
# [11] "perTileQual"                    "perSeqQual"
# [13] "perBaseContent"                 "GCcontent"
# [15] "Ncontent"                       "SeqLengthDist"
# [17] "SeqDuplication"                 "OverrepSeqs"
# [19] "AdapterContent"                 "KmerContent"
# [21] "SeqLength_R1"                   "percentGC_R1"
# [23] "phred20.21_R1"                  "phred48.49_R1"
# [25] "phred76.77_R1"                  "phred100.101_R1"
# [27] "phredGT30_R1"                   "phredGT35_R1"
# [29] "Adapter50.51_R1"                "Adapter70.71_R1"
# [31] "Adapter88.89_R1"                "SeqLength_R2"
# [33] "percentGC_R2"                   "phred20.21_R2"
# [35] "phred48.49_R2"                  "phred76.77_R2"
# [37] "phred100.101_R2"                "phredGT30_R2"
# [39] "phredGT35_R2"                   "Adapter50.51_R2"
# [41] "Adapter70.71_R2"                "Adapter88.89_R2"
# [43] "ERCCsumLogErr"                  "bamFile"
# [45] "trimmed"                        "numReads"
# [47] "numMapped"                      "numUnmapped"
# [49] "overallMapRate"                 "concordMapRate"
# [51] "totalMapped"                    "mitoMapped"
# [53] "mitoRate"                       "totalAssignedGene"
# [55] "gene_Assigned"                  "gene_Unassigned_Ambiguity"
# [57] "gene_Unassigned_MultiMapping"   "gene_Unassigned_NoFeatures"
# [59] "gene_Unassigned_Unmapped"       "gene_Unassigned_MappingQuality"
# [61] "gene_Unassigned_FragmentLength" "gene_Unassigned_Chimera"
# [63] "gene_Unassigned_Secondary"      "gene_Unassigned_Nonjunction"
# [65] "gene_Unassigned_Duplicate"      "rRNA_rate"
# [67] "Flowcell"                       "hasGenotype"

unique(colData(rse_gene)$PrimaryDx)
# [1] "Schizo"  "Control"

table(colData(rse_gene)$PrimaryDx)
# Control  Schizo
#      34      35

table(colData(rse_gene)$Sex)
#  M
# 69

table(colData(rse_gene)$Race)
# CAUC
#   69

table(colData(rse_gene)$Flowcell)
# HVYTYBBXX HW252BBXX
#        34        35

###############################################################################



############################## Data exploration ###############################

colData(rse_gene)$library_size <- apply(assay(rse_gene), 2, sum)
colData(rse_gene)$log10_library_size <- log10(colData(rse_gene)$library_size)
colData(rse_gene)$detected_num_genes <- apply(assay(rse_gene), 2, function(x) {
    length(x[which(x > 0)])
})

qc_metrics <- c("mitoRate", "rRNA_rate", "overallMapRate", "totalAssignedGene", "concordMapRate", "log10_library_size", "detected_num_genes", "RIN")

## Sample variables of interest
sample_variables <- c("PrimaryDx", "Flowcell")

## Function to create boxplots of QC metrics for groups of samples
QC_boxplots <- function(qc_metric, sample_var) {
    if (sample_var == "PrimaryDx") {
        colors <- c("HVYTYBBXX" = "darkmagenta", "HW252BBXX" = "yellow3")
        violin_width <- 0.7
        jitter_width <- 0.1
        x_label <- "PrimaryDx"
        sample_var_v2 <- "Flowcell"
    } else if (sample_var == "Flowcell") {
        colors <- c("Schizo" = "darkgoldenrod3", "Control" = "turquoise3")
        violin_width <- 0.7
        jitter_width <- 0.1
        x_label <- "Flowcell"
        sample_var_v2 <- "PrimaryDx"
    }

    y_label <- str_replace_all(qc_metric, c("_" = " "))

    data <- data.frame(colData(rse_gene))
    plot <- ggplot(data = data, mapping = aes(x = !!rlang::sym(sample_var), y = !!rlang::sym(qc_metric), color = !!rlang::sym(sample_var_v2))) +
        theme_bw() +
        geom_violin(alpha = 0, size = 0.4, color = "black", width = violin_width) +
        geom_jitter(width = jitter_width, alpha = 0.7, size = 2) +
        geom_boxplot(alpha = 0, size = 0.4, width = 0.1, color = "black") +
        scale_color_manual(values = colors) +
        labs(y = y_label, x = x_label) +
        theme(
            axis.title = element_text(size = (9)),
            axis.text = element_text(size = (8))
        )

    return(plot)
}


## Plotting
for (sample_var in sample_variables) {
    width <- 35
    height <- 30
    i <- 1
    plots <- list()
    for (qc_metric in qc_metrics) {
        plots[[i]] <- QC_boxplots(qc_metric, sample_var)
        i <- i + 1
    }
    plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]], plots[[8]], nrow = 3)
    ggsave(paste(here("plots/10_DEA/QC_boxplots_"), sample_var, ".pdf", sep = ""), width = width, height = height, units = "cm")
}



###############################################################################



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
