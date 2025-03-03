## Feb 2025, Louise Huuki-Myers
## For Reviews round 1, plot broad marker genes in bulk data

library("SummarizedExperiment")
library("here")
library("sessioninfo")
library("ggplot2")
library("tidyverse")
library("ComplexHeatmap")
# library("spatialLIBD")
library("RColorBrewer")

# loading sce object
load(here("processed-data","rse_objects", "rse_gene_Habenula_Pilot.rda"),verbose = TRUE)
#rse_gene

rse_gene$PrimaryDx <- gsub("Schizo", "SCZD", rse_gene$PrimaryDx)

# creating plot directory
plot_dir <- here("plots", "02_bulk_qc", "08_bulk_markers_heatmap")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

# sourcing official color palette
source(file = here("code", "99_paper_figs", "source_colors.R"))
  # bulk_colors and sn_colors

#### Check marker genes ####
# list of marker genes
lit_markers = tibble(
    cellType.target = c("Hb", "Hb","MHb"),
    Symbol = c("POU4F1", "GPR151", "TAC3"),
    source = "Literature",
    gene = rowData(rse_gene)$ensemblID[match(c("POU4F1", "GPR151", "TAC3"), rowData(rse_gene)$Symbol)]
)

#### Mean Ratio Markers ####
load(here("processed-data", "06_deconvolution", "run_Bisque", "marker_stats_top_25_genes.Rdata"), verbose = TRUE)

marker_stats |> count(cellType.target)

markers_MR_top <- marker_stats |>
    filter(rank_ratio <= 5) |>
    select(cellType.target, gene, Symbol, rank_ratio) |>
    mutate(source = "MeanRatio") |>
    bind_rows(lit_markers) |>
    mutate(cellType.target = factor(cellType.target, levels = c("Astrocyte", "Endo","Microglia","Oligo","OPC",
                                                                       "Excit.Thal", "Inhib.Thal","Hb", "LHb","MHb")),
                  in_bulk = gene %in% rowData(rse_gene)$ensemblID) |>
    arrange(cellType.target, rank_ratio)

markers_MR_top |> print(n = 48)

## some genes are missing
markers_MR_top |> filter(Symbol %in% lit_markers$Symbol)

markers_MR_top |> filter(!in_bulk)
# cellType.target gene            Symbol     rank_ratio in_bulk
# <fct>           <chr>           <chr>           <int> <lgl>
# 1 Astrocyte       ENSG00000287704 AC012405.1          3 FALSE
# 2 Microglia       ENSG00000249867 LINC02742           1 FALSE
# 3 Oligo           ENSG00000286279 AC108721.2          2 FALSE
# 4 MHb             ENSG00000070748 CHAT                1 FALSE
# 5 MHb             ENSG00000231671 LINC01307           2 FALSE
# 6 MHb             ENSG00000253236 LINC02143           5 FALSE

markers_MR_top <- markers_MR_top |> filter(in_bulk)

#### Complex Heatmap setup ####
rownames(rse_gene) <- rowData(rse_gene)$ensemblID
bulk_logcounts <- assays(rse_gene)$logcounts[markers_MR_top$gene,]
rownames(bulk_logcounts) <- markers_MR_top$Symbol
colnames(bulk_logcounts) <- rse_gene$BrNum

dim(bulk_logcounts)
# [1] 42 68

marker_z_score <- scale(t(bulk_logcounts))

#### prep complex heatmap annotations  ####

## row sample annotations
sample_anno <- as.data.frame(colData(rse_gene), row.names = NULL) |>
    mutate(sn = BrNum %in% c("Br1092", "Br1204", "Br1469", "Br1735", "Br5555", "Br5558", "Br5639")) |>
    select(tot.Hb, PrimaryDx, sn)

rownames(sample_anno) <- rse_gene$BrNum

## reorder matrix by tot.Hb
sample_anno <- sample_anno |> arrange(tot.Hb)
marker_z_score <- marker_z_score[rownames(sample_anno),]

marker_z_score[1:5, 1:5]

row_ha <- rowAnnotation(df = sample_anno,
                        col = list(PrimaryDx = c(Control = "#2a9d8f", SCZD = "#f77f00"),
                                   sn = c(`TRUE` = "black", `FALSE` = "white")))


## column marker gene annotations

broad_cell_type_colors <- c(sn_colors, "Hb" = "#702963",  "MHb" = "#F33A6A", "LHb" = "#0000FF")[levels(markers_MR_top$cellType.target)]

column_ha <- HeatmapAnnotation(
    Marker_Gene = markers_MR_top$cellType.target,
    source = markers_MR_top$source,
    col = list(Marker_Gene = broad_cell_type_colors,
               source = c(`MeanRatio` = "#FFDA85", `Literature` = "#95E4EE"))
)


# printing
pdf(here(plot_dir, "Hb_bulk_heatmap_MRtop5.pdf"), width = 12, height = 10)
print(Heatmap(marker_z_score,
                      name = "Z Score",
                      col = circlize::colorRamp2(seq(-4, 4,  8/10),
                                                 rev(RColorBrewer::brewer.pal(11, "RdBu"))),
                      cluster_rows = TRUE,
                      cluster_columns = FALSE,
                      right_annotation = row_ha,
                      top_annotation = column_ha,
                      column_split = markers_MR_top$cellType.target
))

dev.off()

pdf(here(plot_dir, "Hb_bulk_heatmap_MRtop5_orderHb.pdf"), width = 12, height = 10)
print(Heatmap(marker_z_score,
                      name = "Z Score",
                      col = circlize::colorRamp2(seq(-4, 4,  8/10),
                                                 rev(RColorBrewer::brewer.pal(11, "RdBu"))),
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      right_annotation = row_ha,
                      top_annotation = column_ha,
                      column_split = markers_MR_top$cellType.target
))

dev.off()


sessioninfo::session_info()
