library("here")
library("SummarizedExperiment")
library("scater")
library("variancePartition")
library("ComplexHeatmap")
library("circlize")
library("ggplot2")
library("sessioninfo")

out_plot <- here("plots", "10_DEA", "03_VariancePartition")



########################### Load rse filtered object ##########################

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

dim(rse_gene)
# [1] 22756    68

dim(colData(rse_gene))
# [1]  68 109

###############################################################################



############################## Variance explained #############################

## Set up qc_metrics and colors for the plot
qc_metrics <- c("PrimaryDx", "AgeDeath", "Flowcell", "mitoRate", "rRNA_rate", "totalAssignedGene", "RIN", "abs_ERCCsumLogErr", "tot.Hb", "tot.Thal", "qSV1", "qSV2", "qSV3", "qSV4", "qSV5", "qSV6", "qSV7", "qSV8")
colors <- c(
    "PrimaryDx" = "#e2ef70", "AgeDeath" = "#8093f1", "Flowcell" = "#ddfff7", "mitoRate" = "turquoise4", "rRNA_rate" = "bisque2", "totalAssignedGene" = "blueviolet", "RIN" = "blue3", "abs_ERCCsumLogErr" = "#06d6a0",
    "tot.Hb" = "palegreen3", "tot.Thal" = "skyblue2",
    "qSV1" = "#800f2f", "qSV2" = "#a4133c", "qSV3" = "#c9184a", "qSV4" = "#ff4d6d", "qSV5" = "#ff758f", "qSV6" = "#FF8FA3", "qSV7" = "#FFB3C1", "qSV8" = "#FFCCD5"
)

message(Sys.time(), " - running getVarianceExplained()")

exp_vars <- getVarianceExplained(rse_gene, variables = qc_metrics, exprs_values = "logcounts")

## Plot density graph for each variable
varience_plot <- plotExplanatoryVariables(exp_vars, nvars_to_plot = Inf)
varience_plot <- varience_plot + scale_colour_manual(values = colors) +
    labs(color = "Variables") +
    theme(element_text(size = 15))

ggsave(
    paste(out_plot, "/", "ExplanatoryVars.pdf", sep = ""),
    varience_plot,
    width = 35,
    height = 25,
    units = "cm"
)

###############################################################################



######################### Correlation between variables #######################

formula <- ~ PrimaryDx + AgeDeath + Flowcell + mitoRate + rRNA_rate + totalAssignedGene + RIN + abs_ERCCsumLogErr +
    qSV1 + qSV2 + qSV3 + qSV4 + qSV5 + qSV6 + qSV7 + qSV8 +
    tot.Hb + tot.Thal

message(Sys.time(), " - running canCorPairs()")

corpairs <- canCorPairs(formula, colData(rse_gene))

pdf(paste(out_plot, "/", "CCA_heatmap.pdf", sep = ""), width = 10, height = 10)
Heatmap(
    corpairs,
    col = colorRamp2(c(0, 1), c("#faf0ca", "#3f37c9"))
)
dev.off()

###############################################################################



############################## Variance partition #############################

## with qSVs
formula <- ~ (1 | PrimaryDx) + AgeDeath + (1 | Flowcell) + mitoRate + rRNA_rate + totalAssignedGene + RIN + abs_ERCCsumLogErr +
    qSV1 + qSV2 + qSV3 + qSV4 + qSV5 + qSV6 + qSV7 + qSV8 +
    tot.Hb + tot.Thal

message(Sys.time(), " - running fitExtractVarPartModel()")

## Loop over each gene to fit model and extract variance explained by each variable
varPart <- fitExtractVarPartModel(assays(rse_gene)$logcounts, formula, colData(rse_gene))
# Warning messages:
# 1: Some predictor variables are on very different scales: consider rescaling

# Sort variables by median fraction of variance explained
vp <- sortCols(varPart)

p <- plotVarPart(vp, col = colors) +
    theme(element_text(size = 15))
ggsave(
    filename = paste(out_plot, "/", "VarPartition.pdf", sep = ""),
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
