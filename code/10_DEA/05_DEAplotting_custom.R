library("here")
library("data.table")
library("dplyr")
library("gplots")
library("scater")
library("EnhancedVolcano")
library("sessioninfo")

out_plot <- here("plots", "10_DEA", "04_DEA")

#### load data ####

gene_de <- read.csv(here("processed-data", "15_compare_vs_BrainSeq", "TableSxx_SZCD_vs_Control_gene_results.csv"))
head(gene_de)

FDR_cut = 10e-02

keyvals <- ifelse(
   gene_de$Hb_adj.P.Val >= FDR_cut, "#6c757d", ifelse(
       gene_de$Hb_logFC < 0, "#2a9d8f", "#f77f00"
    )
)
names(keyvals)[keyvals == "#2a9d8f"] <- paste0("Down - FDR < ", FDR_cut)
names(keyvals)[keyvals == "#6c757d"] <- "Not significant"
names(keyvals)[keyvals == "#f77f00"] <- paste0("Up - FDR < ", FDR_cut)

hval <- max((gene_de |> dplyr::filter(Hb_adj.P.Val < 0.1 & Hb_adj.P.Val >= 0.09)) |> pull(Hb_P.Value))
hval <- format(hval, scientific = TRUE)

y_lim <- ceiling(max(-log10(gene_de$Hb_P.Value)))

select_genes <- c("CCDC141", "QPRT", "HES5", "EHMT2", "NOTCH4", "DRD5", "GFRA1", "NDST1" ,"FLRT1")

gene_volcano <- EnhancedVolcano(gene_de,
                x = "Hb_logFC",
                y = "Hb_P.Value",
                selectLab = select_genes,
                labSize = 6.0,
                drawConnectors = TRUE,
                arrowheads = FALSE,
                max.overlaps = Inf,
                labCol = "black",
                pCutoff = hval,
                FCcutoff = 1,
                lab = gene_de$Symbol,
                colCustom = keyvals,
                caption = paste0(nrow(gene_de), " total genes\n",
                                 sum(gene_de$Hb_adj.P.Val < FDR_cut, na.rm = TRUE),
                                 " significant genes"),
                title = NULL,
                subtitle = NULL) +
    ylim(c(0, y_lim)) +
    xlim(c(-3, 3)) +
    ylab("-log10(p-value)") +
    xlab("log2FC (SCZD vs Control)")

ggsave(gene_volcano, filename = here(out_plot, "VolcanoPlot-gene_custom.png"), height = 9)
ggsave(gene_volcano, filename = here(out_plot, "VolcanoPlot-gene_custom.pdf"), height = 9)

######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################

