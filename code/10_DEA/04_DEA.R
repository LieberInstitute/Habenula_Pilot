library("here")
library("SummarizedExperiment")
library("edgeR")
library("limma")
library("dplyr")
library("gplots")
library("scater")
library("EnhancedVolcano")
library("sessioninfo")

out_plot <- here("plots", "10_DEA", "04_DEA")
out_data <- here("processed-data", "10_DEA", "04_DEA")



############################## Load rse objects ###############################

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



######################## Functions for DEA and plotting #######################

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

    all_df <- top_genes
    all_df$ensemblID <- NULL
    all_df <- tibble::rownames_to_column(all_df, "ensemblID")
    write.table(all_df, file = paste0(out_data, "/DEA_AllGenes_", model_name, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

    return(top_genes)
}

## Function to make volcano plot
plot_volc <- function(top_genes, FDR_cut, model_name, hval) {
    sig_df <- top_genes %>% filter(adj.P.Val < FDR_cut)
    sig_df$ensemblID <- NULL
    sig_df <- tibble::rownames_to_column(sig_df, "ensemblID")

    write.table(sig_df %>% filter(adj.P.Val < FDR_cut), file = paste0(out_data, "/DEA_SigGenes_FDR", gsub(as.character(FDR_cut), pattern = "0\\.", replacement = ""), "_", model_name, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

    ## Format data
    outGenes_plot <- top_genes %>%
        select(logFC, P.Value, adj.P.Val, ensemblID, Symbol)
    rownames(outGenes_plot) <- uniquifyFeatureNames(
        outGenes_plot$ensemblID,
        outGenes_plot$Symbol
    )

    ## Select colors
    keyvals <- ifelse(
        outGenes_plot$adj.P.Val >= FDR_cut, "#6c757d", ifelse(
            outGenes_plot$logFC < 0, "#2a9d8f", "#f77f00")
    )
    names(keyvals)[keyvals == "#2a9d8f"] <- paste0("Down - FDR < ", FDR_cut)
    names(keyvals)[keyvals == "#6c757d"] <- "Not significant"
    names(keyvals)[keyvals == "#f77f00"] <- paste0("Up - FDR < ", FDR_cut)

    ## Genes to highlight
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
        caption = paste0(nrow(outGenes_plot), " total genes\n", sum(outGenes_plot$adj.P.Val < FDR_cut), " significant genes"),
        title = NULL,
        subtitle = NULL
    ) + ylim(c(0, 6)) +
        xlim(c(-3, 3)) +
        ylab("-log10(p-value)") +
        xlab("log2FC (SCZD vs Control)")

    pdf(paste0(out_plot, "/VolcanoPlot_", model_name, ".pdf"),
        height = 10,
        width = 8
    )
    print(volcano_plot)
    dev.off()
}

###############################################################################



##################################### DEA #####################################

formula <- ~ PrimaryDx + AgeDeath + Flowcell + mitoRate + rRNA_rate + RIN + totalAssignedGene + abs_ERCCsumLogErr +
    qSV1 + qSV2 + qSV3 + qSV4 + qSV5 + qSV6 + qSV7 + qSV8 +
    tot.Hb + tot.Thal

res_formula <- DE_analysis(
    rse_gene = rse_gene,
    formula = formula,
    coef = "PrimaryDxSchizo",
    model_name = "qc-totAGene-qSVs-Hb-Thal"
)

hval <- max((res_formula %>% dplyr::filter(adj.P.Val < 0.1 & adj.P.Val >= 0.09))$P.Value)
hval <- format(hval, scientific = TRUE)

plot_volc(res_formula, FDR_cut = 10e-02, model_name = "qc-totAGene-qSVs-Hb-Thal", hval = as.numeric(hval))

###############################################################################



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
