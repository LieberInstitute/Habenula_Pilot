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

DE_all_files <- list.files(
    here(
        "processed-data",
        "10_DEA",
        "04_DEA"
    ),
    pattern = "DEA_All*",
    full.names = TRUE
)
DE_all_files <- DE_all_files[1:3]

DE_all <- lapply(DE_all_files, fread, data.table = FALSE)

DE_sig_files <- list.files(
    here(
        "processed-data",
        "10_DEA",
        "04_DEA"
    ),
    pattern = "DEA_Sig*",
    full.names = TRUE
)
DE_sig_files <- DE_sig_files[1:3]

DE_sig <- lapply(DE_sig_files, fread, data.table = FALSE)

###############################################################################



######################### Prepare data fo GO analysis #########################

## Converting ensembl IDs to Entrez IDs in case some of my data.frames does not include that column
DE_all <- lapply(DE_all, function(x) {
    entrez <- mapIds(org.Hs.eg.db, keys = x$ensemblID, keytype = "ENSEMBL", column = "ENTREZID")
    entrez[is.na(names(entrez))] <- NA
    entrez <- unlist(entrez)
    x$EntrezID <- entrez
    return(x)
})

DE_sig <- lapply(DE_sig, function(x) {
    entrez <- mapIds(org.Hs.eg.db, keys = x$ensemblID, keytype = "ENSEMBL", column = "ENTREZID")
    entrez[is.na(names(entrez))] <- NA
    entrez <- unlist(entrez)
    x$EntrezID <- entrez
    return(x)
})


###############################################################################



##################### Run GO and KEGG enrichment analysis #####################

go <- compareCluster(sigGene,
    fun = "enrichGO",
    universe = geneUniverse,
    OrgDb = org.Hs.eg.db,
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.15,
    #qvalueCutoff = 0.05,
    readable = TRUE
)

kegg <- compareCluster(sigGene,
    fun = "enrichKEGG",
    universe = geneUniverse,
    organism = "human",
    pvalueCutoff = 0.15,
    qvalueCutoff = 0.05
)

###############################################################################



########################### Plot GO and KEGG results ##########################

## Function to plot GO results
plot_go <- function(ont, title_p, path, filename, size) {
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

    ggsave(filename = filename, path = path, dotplot_1, height = size[1], width = size[2])
}

plot_go(ont = "BP", title_p = "Biological Process", filename = "GOenrichment_BP_qc-snpPCs-Hb.pdf", path = out_plot, size = c(4, 8))
plot_go(ont = "CC", title_p = "Cellular Component", filename = "GOenrichment_CC_qc-snpPCs-Hb.pdf", path = out_plot, size = c(6, 5))
plot_go(ont = "MF", title_p = "Molecular Function", filename = "GOenrichment_MF_qc-snpPCs-Hb.pdf", path = out_plot, size = c(6, 6.5))

## Plot KEGG results
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

ggsave(filename = "KEGGenrichment_qc-snpPCs-Hb.pdf", path = out_plot, dotplot_1, height = 6, width = 5)

###############################################################################



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
