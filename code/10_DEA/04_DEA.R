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

rse_objects <- list.files(
    here(
        "processed-data",
        "rse_objects"
    ),
    full.names = TRUE
)

lapply(rse_objects, load, verbose = TRUE, .GlobalEnv)
# Loading objects:
#   rse_exon
# Loading objects:
#   rse_gene
# Loading objects:
#   rse_jx
# Loading objects:
#   rse_tx

###############################################################################



######################## Functions for DEA and plotting #######################

## Function to do all the DEA analysis with limma
DE_analysis <- function(rse, formula, coef, model_name, FDR_cut = 10e-02, rse_type) {
    pdf(file = paste0(out_plot, "/summary-", rse_type, "_", model_name, ".pdf"))
    # par(mfrow = c(3, 2))
    par(cex = 0.7, mai = c(0.1, 0.1, 0.1, 0.1))

    ## Model matrix
    model <- model.matrix(formula, data = colData(rse))

    if (rse_type != "tx") {
        ## Use previous norm factors to scale the raw library sizes
        rse_scaled <- calcNormFactors(rse)

        par(fig = c(0.05, 0.5, 0.55, 0.95))
        ## Transform counts to log2(CPM): estimate mean-variance relationship for
        ## each gene
        vFeat <- voom(rse_scaled, design = model, plot = TRUE)

        ## Fit linear model for each gene
        fitFeat <- lmFit(vFeat)
    } else {
        fitFeat <- lmFit(assays(rse)$logcounts, design = model)
    }

    ## Empirical Bayesian calculation to obtain our significant genes: compute
    ## moderated F and t-statistics, and log-odds of DE
    eBFeat <- eBayes(fitFeat)

    par(fig = c(0.55, 1, 0.55, 0.95), new = TRUE)
    ## Plot average log expression vs logFC
    limma::plotMA(eBFeat,
        coef = coef, xlab = "Mean of normalized counts",
        ylab = "logFC"
    )

    par(fig = c(0.05, 0.5, 0.1, 0.5), new = TRUE)

    ## Plot -log(p-value) vs logFC
    volcanoplot(eBFeat, coef = coef)

    ## Select top-ranked genes
    top_genes <- topTable(eBFeat, coef = coef, p.value = 1, number = nrow(rse), sort.by = "none")

    ## Histogram of adjusted p values
    par(fig = c(0.55, 1, 0.1, 0.5), new = TRUE)
    hist(top_genes$P.Value, xlab = "p-value", main = "")

    par(fig = c(0.25, 0.75, 0, 0.1), new = TRUE)
    textplot(capture.output(summary(top_genes$adj.P.Val)))

    par(fig = c(0.0, 0.25, 0, 0.1), new = TRUE)
    textplot(capture.output(length(unique(top_genes$adj.P.Val))))

    dev.off()

    df_class <- sapply(top_genes, class) == "list"
    if (sum(df_class) > 0) {
        for (i in which(df_class)) {
            top_genes[, i] <- unlist(lapply(top_genes[, i], paste, collapse = ","))
        }
    }

    ## Save all genes/tx/jx/exons and then only the significant ones
    all_df <- top_genes
    all_df <- tibble::rownames_to_column(all_df, "r.name")

    write.table(all_df, file = paste0(out_data, "/DEA_All-", rse_type, "_", model_name, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

    sig_df <- top_genes %>% filter(adj.P.Val < FDR_cut)
    sig_df <- tibble::rownames_to_column(sig_df, "r.name")

    write.table(sig_df %>% filter(adj.P.Val < FDR_cut), file = paste0(out_data, "/DEA_Sig-", rse_type, "_FDR", gsub(as.character(FDR_cut), pattern = "0\\.", replacement = ""), "_", model_name, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

    return(top_genes)
}

###############################################################################



##################################### DEA #####################################

formula <- ~ PrimaryDx + AgeDeath + Flowcell + mitoRate + rRNA_rate + RIN + totalAssignedGene + abs_ERCCsumLogErr +
    qSV1 + qSV2 + qSV3 + qSV4 + qSV5 + qSV6 + qSV7 + qSV8 +
    tot.Hb + tot.Thal

## rse_gene
res_formula <- DE_analysis(
    rse = rse_gene,
    formula = formula,
    coef = "PrimaryDxSchizo",
    model_name = "qc-totAGene-qSVs-Hb-Thal",
    rse_type = "gene"
)

## rse_jx
res_formula <- DE_analysis(
    rse = rse_jx,
    formula = formula,
    coef = "PrimaryDxSchizo",
    model_name = "qc-totAGene-qSVs-Hb-Thal",
    rse_type = "jx"
)
# Warning message:
# Partial NA coefficients for 8 probe(s)

## rse_exon
res_formula <- DE_analysis(
    rse = rse_exon,
    formula = formula,
    coef = "PrimaryDxSchizo",
    model_name = "qc-totAGene-qSVs-Hb-Thal",
    rse_type = "exon"
)

## rse_tx
res_formula <- DE_analysis(
    rse = rse_tx,
    formula = formula,
    coef = "PrimaryDxSchizo",
    model_name = "qc-totAGene-qSVs-Hb-Thal",
    rse_type = "tx"
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
