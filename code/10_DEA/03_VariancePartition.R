library("here")
library("SummarizedExperiment")
library("scater")
library("variancePartition")
library("ComplexHeatmap")
library("circlize")
library("ggplot2")
library("sessioninfo")

output_path <- here("plots", "10_DEA", "03_VariancePartition")



########################### Load rse filtered object ##########################

load(
    here(
        "processed-data",
        "10_DEA",
        "rse_gene_filt.Rdata"
    ),
    verbose = TRUE
)
# Loading objects:
#   rse_gene_filt

lobstr::obj_size(rse_gene_filt)
# 27.27 MB

class(rse_gene_filt)
# [1] "RangedSummarizedExperiment"
# attr(,"package")
# [1] "SummarizedExperiment"

###############################################################################



############################## Variance explained #############################

## Set up qc_metrics and colors for the plot
qc_metrics <- c("mitoRate", "rRNA_rate", "overallMapRate", "totalAssignedGene", "concordMapRate", "library_size", "detected_num_genes", "RIN", "abs_ERCCsumLogErr", "PrimaryDx", "Flowcell", "AgeDeath")
colors <- c("mitoRate" = "turquoise4", "rRNA_rate" = "bisque2", "overallMapRate" = "indianred1", "totalAssignedGene" = "blueviolet", "concordMapRate" = "lightsalmon", "library_size" = "palegreen3", "detected_num_genes" = "skyblue2", "RIN" = "blue3", "abs_ERCCsumLogErr" = "#06d6a0", "PrimaryDx" = "#a14a76", "Flowcell" = "#fdc500", "AgeDeath" = "#dda15e")

exp_vars <- getVarianceExplained(rse_gene_filt, variables = qc_metrics, exprs_values = "logcounts")

## Plot density graph for each variable
varience_plot <- plotExplanatoryVariables(exp_vars, theme_size = 12, nvars_to_plot = Inf)
varience_plot <- varience_plot + scale_colour_manual(values = colors) +
    labs(color = "Variables")

ggsave(
    paste(output_path, "/", "ExplanatoryVars.pdf", sep = ""),
    varience_plot,
    width = 35,
    height = 25,
    units = "cm"
)

###############################################################################



######################### Correlation between variables #######################

formula <- ~ PrimaryDx + AgeDeath + Flowcell + mitoRate + rRNA_rate + overallMapRate + totalAssignedGene + concordMapRate + library_size + detected_num_genes + RIN + abs_ERCCsumLogErr

corpairs <- canCorPairs(formula, colData(rse_gene_filt))

pheatmap(
    corpairs,
    color = hcl.colors(50, "YlOrRd", rev = TRUE),
    fontsize = 8,
    border_color = "black",
    height = 6,
    width = 6.5,
    filename = paste(output_path, "/", "CCA_heatmap.pdf", sep = "")
)

###############################################################################



############################## Variance partition #############################

formula <- ~ (1 | PrimaryDx) + AgeDeath + (1 | Flowcell) + mitoRate + rRNA_rate + overallMapRate + totalAssignedGene + concordMapRate + library_size + detected_num_genes + RIN + abs_ERCCsumLogErr

## Loop over each gene to fit model and extract variance explained by each variable
varPart <- fitExtractVarPartModel(assays(rse_gene_filt)$logcounts, formula, colData(rse_gene_filt))
# Warning messages:
# 1: Some predictor variables are on very different scales: consider rescaling

# Sort variables by median fraction of variance explained
vp <- sortCols(varPart)

p <- plotVarPart(vp)
ggsave(
    filename = paste(output_path, "/", "VarPartition.pdf", sep = ""),
    p,
    width = 40,
    height = 20,
    units = "cm"
)

###############################################################################



######################## Variance partition - filtered ########################

## Plots without overallMapRate, concordMapRate, detected_num_genes
formula <- ~ (1 | PrimaryDx) + AgeDeath + (1 | Flowcell) + mitoRate + rRNA_rate + totalAssignedGene + RIN + abs_ERCCsumLogErr + library_size

varPart <- fitExtractVarPartModel(assays(rse_gene_filt)$logcounts, formula, colData(rse_gene_filt))
# Warning messages:
# 1: Some predictor variables are on very different scales: consider rescaling

# Sort variables by median fraction of variance explained
vp <- sortCols(varPart)

p <- plotVarPart(vp)
ggsave(
    filename = paste(output_path, "/", "VarPartition_filtered.pdf", sep = ""),
    p,
    width = 40,
    height = 20,
    units = "cm"
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
