library("here")
library("data.table")
library("dplyr")
library("clusterProfiler")
library("org.Hs.eg.db")
library("ggplot2")
library("sessioninfo")

out_plot <- here("plots", "10_DEA", "07_GOplot")
if (!dir.exists(out_plot)) dir.create(out_plot)



########################### Load GO and KEGG results ##########################

go_results_files <- list.files(
    here(
        "processed-data",
        "10_DEA",
        "05_GOenrichment"
    ),
    pattern = "go_.+.csv",
    full.names = TRUE
)
go_results <- lapply(go_results_files, fread, data.table = FALSE)
names(go_results) <- c("exon", "gene", "jx")

kegg_results_files <- list.files(
    here(
        "processed-data",
        "10_DEA",
        "05_GOenrichment"
    ),
    pattern = "kegg_.+.csv",
    full.names = TRUE
)
kegg_results <- lapply(kegg_results_files, fread, data.table = FALSE)
names(kegg_results) <- c("exon", "gene", "jx")

###############################################################################



################################ Filter results ###############################

go_results_filt <- lapply(go_results, function(feature) {
    feature_filt <- feature %>%
        filter(qvalue < 0.2) %>%
        arrange(qvalue)
    return(feature_filt)
})
lapply(go_results_filt, dim)
# $exon
# [1] 115  12
# $gene
# [1] 36 12
# $jx
# [1] 1617   12

kegg_results_filt <- lapply(kegg_results, function(feature) {
    feature_filt <- feature %>%
        filter(qvalue < 0.2) %>%
        arrange(qvalue)
    return(feature_filt)
})
lapply(kegg_results_filt, dim)
# $exon
# [1]  1 11
# $gene
# [1]  0 11
# $jx
# [1] 96 11

###############################################################################



########################### Plot GO and KEGG results ##########################
go = go_results_filt$exon
ont = "BP"

## Function to plot GO results
plot_go <- function(go, ont, title_p, path, filename, size) {
    leng_ont <- dim(filter(go, ONTOLOGY == ont))[1]
     ifelse(leng_ont >= 10,
         go_mat <- filter(go, ONTOLOGY == ont)[1:10,],
         go_mat <- filter(go, ONTOLOGY == ont)[1:leng_ont,]
     )
    #go_mat <- filter(go, ONTOLOGY == ont)[1:leng_ont,]
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

## Exon
plot_go(go = go_results_filt$exon, ont = "BP", title_p = "Biological Process", filename = "GO-BP_exon.pdf", path = out_plot, size = c(4, 6))
plot_go(go = go_results_filt$exon, ont = "CC", title_p = "Cellular Component", filename = "GO-CC_exon.pdf", path = out_plot, size = c(6, 5))
plot_go(go = go_results_filt$exon, ont = "MF", title_p = "Molecular Function", filename = "GO-MF_exon.pdf", path = out_plot, size = c(6, 10))

## Gene
plot_go(go = go_results_filt$gene, ont = "BP", title_p = "Biological Process", filename = "GO-BP_gene.pdf", path = out_plot, size = c(4, 6))
plot_go(go = go_results_filt$gene, ont = "CC", title_p = "Cellular Component", filename = "GO-CC_gene.pdf", path = out_plot, size = c(5, 8))
# plot_go(go = go_results_filt$gene, ont = "MF", title_p = "Molecular Function", filename = "GO-MF_gene.pdf", path = out_plot, size = c(6, 10))

## Jx
plot_go(go = go_results_filt$jx, ont = "BP", title_p = "Biological Process", filename = "GO-BP_jx.pdf", path = out_plot, size = c(4, 6))
plot_go(go = go_results_filt$jx, ont = "CC", title_p = "Cellular Component", filename = "GO-CC_jx.pdf", path = out_plot, size = c(6, 5))
plot_go(go = go_results_filt$jx, ont = "MF", title_p = "Molecular Function", filename = "GO-MF_jx.pdf", path = out_plot, size = c(6, 10))


## Plot KEGG results
dotplot_1 <- ggplot(kegg_results_filt$jx, aes(Cluster, Description)) +
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

ggsave(filename = "KEGG_jx.pdf", path = out_plot, dotplot_1, height = 10, width = 8)

###############################################################################



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################

