library("here")
library("data.table")
library("dplyr")
library("clusterProfiler")
library("org.Hs.eg.db")
library("ggplot2")
library("sessioninfo")

out_plot <- here("plots", "10_DEA", "05_GOenrichment")
out_data <- here("processed-data", "10_DEA", "05_GOenrichment")



####################### Load tsv with all genes from DEA ######################

## For now I'm just gonna load these two models: qc-snpPCs-Hb and
## qc-totAssGene-snpPCs-Hb beacuse they had the best results

DE_qc_snpPCs_Hb <- fread(
    here(
        "processed-data",
        "10_DEA",
        "04_DEA",
        "DEA_AllGenes_qc-snpPCs-Hb.tsv"
    ),
    sep = "\t",
    data.table = FALSE,
    stringsAsFactors = FALSE
)

DE_qc_totAssGene_snpPCs_Hb <- fread(
    here(
        "processed-data",
        "10_DEA",
        "04_DEA",
        "DEA_AllGenes_qc-totAssGene-snpPCs-Hb.tsv"
    ),
    sep = "\t",
    data.table = FALSE,
    stringsAsFactors = FALSE
)

###############################################################################


######################### Prepare data fo GO analysis #########################

sigGene <- DE_qc_snpPCs_Hb %>% filter(adj.P.Val < 0.15 & abs(logFC) > 1)

sigGene <- split(sigGene$EntrezID, sign(sigGene$logFC))
sigGene <- lapply(sigGene, function(x) x[!is.na(x)])

geneUniverse <- as.character(DE_qc_snpPCs_Hb$EntrezID)
geneUniverse <- geneUniverse[!is.na(geneUniverse)]

###############################################################################



######################### Run GO and KEGG enrichment analysis #########################

go <- compareCluster(sigGene,
    fun = "enrichGO",
    universe = geneUniverse,
    OrgDb = org.Hs.eg.db,
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.05,
    readable = TRUE
)

kegg <- compareCluster(sigGene,
    fun = "enrichKEGG",
    universe = geneUniverse,
    organism = "human",
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.05
)


plot_go <- function(ont, title_p, path, filename) {
    leng_ont <- dim(filter(go, ONTOLOGY == ont))[1]
    ifelse(leng_ont >= 10,
        go_mat <- filter(go, ONTOLOGY == ont)[1:10],
        go_mat <- filter(go, ONTOLOGY == ont)[1:leng_ont]
        )
    dotplot_1 <- ggplot(go_mat, aes(Cluster, Description)) +
        theme_bw() +
        geom_point(aes(color = p.adjust, size = Count)) +
        scale_color_gradientn(
            colours = c("#f7ca64", "#46bac2", "#7e62a3"),
            trans = "log10",
            guide = guide_colorbar(reverse = TRUE, order = 1)
        ) +
        scale_size_continuous(range = c(2, 10)) +
        xlab("Cluster") +
        ylab("") +
        ggtitle(title_p)

    ggsave(filename = filename, path = path, dotplot_1, height = 6, width = 10)
}

plot_go(ont = "BP", title_p = "Biological Process", filename = "GOenrichment_BP.pdf", path = out_plot)

plot_go(ont = "CC", title_p = "Cellular Component", filename = "GOenrichment_CC.pdf", path = out_plot)

plot_go(ont = "MF", title_p = "Molecular Function", filename = "GOenrichment_MF.pdf", path = out_plot)



dotplot_1 <- ggplot(kegg, aes(Cluster, Description)) +
    theme_bw() +
    geom_point(aes(color = p.adjust, size = Count)) +
    scale_color_gradientn(
        colours = c("#f7ca64", "#46bac2", "#7e62a3"),
        trans = "log10",
        guide = guide_colorbar(reverse = TRUE, order = 1)
    ) +
    scale_size_continuous(range = c(2, 10)) +
    xlab("Cluster") +
    ylab("") +
    ggtitle("KEGG")

ggsave(filename = "KEGGenrichment.pdf", path = out_plot, dotplot_1, height = 6, width = 5)


###############################################################################



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
