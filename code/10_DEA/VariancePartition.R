library("here")
library("SummarizedExperiment")
library("ggplot2")
library("scater")
library("rlang")
library("cowplot")
library("Hmisc")
library("lme4")
library("variancePartition")
library("reshape2")
library("pheatmap")
library("sessioninfo")



########################### Load rse filtered object ##########################

load(
    here(
        "processed-data", "10_DEA", "rse_gene_filt.Rdata"
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
qc_metrics <- c("mitoRate", "rRNA_rate", "overallMapRate", "totalAssignedGene", "concordMapRate", "log10_library_size", "detected_num_genes", "RIN")
colors=c("mitoRate"="turquoise4", "rRNA_rate" = "bisque2", "overallMapRate"="indianred1", "totalAssignedGene" = "blueviolet", "concordMapRate"="lightsalmon", "log10_library_size"="palegreen3", "detected_num_genes"="skyblue2", "RIN"="blue3")

exp_vars<-getVarianceExplained(rse_gene_filt, variables=qc_metrics, exprs_values = "logcounts")

## Plot density graph for each variable
varience_plot <- plotExplanatoryVariables(exp_vars, theme_size = 12, nvars_to_plot = Inf)
varience_plot <- varience_plot + scale_colour_manual(values = colors) +
    labs(color="Variables")

ggsave(filename = here("plots/10_DEA/ExplanatoryVars.pdf"), varience_plot, width = 35, height = 25, units = "cm")

###############################################################################


######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################

