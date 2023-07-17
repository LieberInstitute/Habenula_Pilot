library("here")
library("SummarizedExperiment")
library("edgeR")
library("limma")
library("dplyr")
library("gplots")
library("EnhancedVolcano")
library("sessioninfo")

out_plot <- here("plots", "10_DEA", "04_DEA")
out_data <- here("processed-data", "10_DEA", "04_DEA")



###################### Load rse and rse filtered objects ######################

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

###############################################################################



#################### E### Functions for DEA and plotting #######################


## Function to do all the DEA analysis with limma
DE_analysis <- function(rse_gene, formula, coef, model_name) {
    pdf(file = paste0(out_plot, "/DEAplots_", model_name, ".pdf"))
    # par(mfrow = c(3, 2))
    par(cex = 0.7, mai = c(0.1, 0.1, 0.1, 0.1))

    ## Model matrix
    model <- model.matrix(formula, data = colData(rse_gene))

    ## Use previous norm factors to scale the raw library sizes
    RSE_scaled <- calcNormFactors(rse_gene)

    par(fig = c(0.05, 0.5, 0.55, 0.95))
    ## Transform counts to log2(CPM): estimate mean-variance relationship for
    ## each gene
    vGene <- voom(RSE_scaled, design = model, plot = TRUE)

    ## Fit linear model for each gene
    fitGene <- lmFit(vGene)

    ## Empirical Bayesian calculation to obtain our significant genes: compute
    ## moderated F and t-statistics, and log-odds of DE
    eBGene <- eBayes(fitGene)

    par(fig = c(0.55, 1, 0.55, 0.95), new = TRUE)
    ## Plot average log expression vs logFC
    limma::plotMA(eBGene,
        coef = coef, xlab = "Mean of normalized counts",
        ylab = "logFC"
    )

    par(fig = c(0.05, 0.5, 0.1, 0.5), new = TRUE)

    ## Plot -log(p-value) vs logFC
    volcanoplot(eBGene, coef = coef)

    ## Select top-ranked genes
    top_genes <- topTable(eBGene, coef = coef, p.value = 1, number = nrow(rse_gene), sort.by = "none")

    ## Histogram of adjusted p values
    par(fig = c(0.55, 1, 0.1, 0.5), new = TRUE)
    hist(top_genes$P.Value, xlab = "p-value", main = "")

    par(fig = c(0.25, 0.75, 0, 0.1), new = TRUE)
    textplot(capture.output(summary(top_genes$adj.P.Val)))

    par(fig = c(0.0, 0.25, 0, 0.1), new = TRUE)
    textplot(capture.output(length(unique(top_genes$adj.P.Val))))

    dev.off()

    write.table(top_genes, file = paste0(out_data, "/DEA_AllGenes_", model_name, ".tsv"), sep = "\t", quote = FALSE)

    return(top_genes)
}

## Function to make volcano plot
plot_volc <- function(top_genes, FDR_cut, model_name, hval) {
    outGenes_plot <- top_genes %>%
        select(logFC, P.Value, adj.P.Val, ensemblID, Symbol)

    keyvals <- ifelse(
        outGenes_plot$adj.P.Val > FDR_cut, "#f2e8cf", "#a7c957"
    )

    names(keyvals)[keyvals == "#a7c957"] <- paste0("FDR < ", FDR_cut)
    names(keyvals)[keyvals == "#f2e8cf"] <- "Not significant"

    select <- outGenes_plot %>% filter(adj.P.Val < FDR_cut & abs(logFC) > 1)

    volcano_plot <- EnhancedVolcano(outGenes_plot,
        x = "logFC",
        y = "P.Value",
        selectLab = select$Symbol,
        labSize = 6.0,
        drawConnectors = TRUE,
        arrowheads = FALSE,
        max.overlaps = Inf,
        labCol = "black",
        pCutoff = hval,
        FCcutoff = 1,
        lab = rownames(outGenes_plot),
        colCustom = keyvals,
        caption = paste0("total = ", nrow(outGenes_plot), " genes"),
        title = "",
        subtitle = ""
    ) + ylim(c(0, 6)) +
        xlim(c(-3, 3))

    pdf(paste0(out_plot, "/VolcanoPlot_", model_name, ".pdf"),
        height = 10,
        width = 8
    )
    print(volcano_plot)
    dev.off()
}

###############################################################################



##################################### DEA #####################################


plot_volc(res_formula,
    FDR_cut = 10e-02,
    model_name = "qc-totAssGene-snpPCs-Hb-Thal-qSVs"
)

summary(res_formula$adj.P.Val)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.07587 0.79389 0.90266 0.86291 0.96492 0.99996

###############################################################################



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
