library("here")
library("data.table")
library("dplyr")
library("clusterProfiler")
library("org.Hs.eg.db")
library("ggplot2")
library("sessioninfo")

out_plot <- here("plots", "10_DEA", "07_GOplot")
if (!dir.exists(out_plot)) dir.create(out_plot)



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

