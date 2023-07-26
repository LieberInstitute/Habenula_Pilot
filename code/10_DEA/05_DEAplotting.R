library("here")
library("data.table")
library("dplyr")
library("gplots")
library("scater")
library("EnhancedVolcano")
library("sessioninfo")

out_plot <- here("plots", "10_DEA", "04_DEA")



################## Load significant results from DEA objects ##################

dea_res <- list.files(
    here(
        "processed-data", "10_DEA", "04_DEA"
    ),
    pattern = "DEA_All*",
    full.names = TRUE
)

dea_res <- lapply(dea_res, fread, data.table = FALSE)


###############################################################################



######################### Function to make volcanoplot ########################

plot_volc <- function(dea_feat_res, FDR_cut = 10e-02, model_name, hval, rse_type) {
    ## Format data
    outFeat_plot <- dea_feat_res %>%
    select(logFC, P.Value, adj.P.Val, r.name, Symbol)
rownames(outFeat_plot) <- uniquifyFeatureNames(
    outFeat_plot$r.name,
    outFeat_plot$Symbol
)

    ## Select colors
    keyvals <- ifelse(
        outFeat_plot$adj.P.Val >= FDR_cut, "#6c757d", ifelse(
            outFeat_plot$logFC < 0, "#2a9d8f", "#f77f00"
        )
    )
    names(keyvals)[keyvals == "#2a9d8f"] <- paste0("Down - FDR < ", FDR_cut)
    names(keyvals)[keyvals == "#6c757d"] <- "Not significant"
    names(keyvals)[keyvals == "#f77f00"] <- paste0("Up - FDR < ", FDR_cut)

    ## Set up caption
    if (rse_type == "gene") {
        capt_txt <- paste0(nrow(outFeat_plot), " total genes\n", sum(outFeat_plot$adj.P.Val < FDR_cut, na.rm = TRUE), " significant genes")
    }
    if (rse_type == "jx") {
        capt_txt <- paste0(nrow(outFeat_plot), " total jx\n", sum(outFeat_plot$adj.P.Val < FDR_cut, na.rm = TRUE), " significant jx")
    }
    if (rse_type == "exon") {
        capt_txt <- paste0(nrow(outFeat_plot), " total exons\n", sum(outFeat_plot$adj.P.Val < FDR_cut, na.rm = TRUE), " significant exons")
    }

    y_lim <- ceiling(max(-log10(outFeat_plot$P.Value)))

    ## Features to highlight
    select <- outFeat_plot %>% filter(adj.P.Val < FDR_cut & abs(logFC) > 1)

    volcano_plot <- EnhancedVolcano(outFeat_plot,
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
        lab = rownames(outFeat_plot),
        colCustom = keyvals,
        caption = capt_txt,
        title = NULL,
        subtitle = NULL
    ) + ylim(c(0, y_lim)) +
        xlim(c(-3, 3)) +
        ylab("-log10(p-value)") +
        xlab("log2FC (SCZD vs Control)")

    png(paste0(out_plot, "/VolcanoPlot-", rse_type, "_", model_name, ".png"),
        height = 10*96,
        width = 8*96
    )
    print(volcano_plot)
    dev.off()
}

###############################################################################



############################## Plot volcano plots #############################

mapply(function(dea_feat_res, rse_type) {
    hval <- max((dea_feat_res %>% dplyr::filter(adj.P.Val < 0.1 & adj.P.Val >= 0.09))$P.Value)
    hval <- format(hval, scientific = TRUE)

    plot_volc(dea_feat_res, model_name = "qc-totAGene-qSVs-Hb-Thal", hval = as.numeric(hval), rse_type = rse_type)
}, dea_res, c("exon", "gene", "jx", "tx"))

###############################################################################



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################

