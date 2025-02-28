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
official_markers = list(
  "Oligo" = c("MOBP"),
  "OPC" = c("PDGFRA"),
  "Micro" = c("CSF1R"),
  "Astro" = c("AQP4"),
  "Endo" = c("ITIH5"),
  "Thal" = c("LYPD6B", "ADARB2", "RORB"),
  # "LHb.1" = c("LINC02653"), #  , ATP8B1
  # "LHb.2" = c("AC073071.1"),
  # "LHb.3" = c ("ENTHD1"),
  # "LHb.4" = c("TLE2"),
  # "LHb.5" = c("LINC01619"),
  # "LHb.6" = c("TACR3"),
  # "LHb.7" = c("AC008619.1"),
  # "MHb.1" = c("EXOC1L"),
  # "MHb.2" = c("CHAT"),
  # "MHb.3" = c("BHLHE22"),
  'Hb' = c("POU4F1", "GPR151", "TAC3"), # BARHL1
  "MHb" = c("CHRNB4"),
  "LHb" = c("HTR2C"),
  'Neuron' = c('SYT1'),
  'Exc_Neuron' = c('SLC17A6'),
  'Inh_Neuron' = c('GAD1')
)

marker_stats <- readxl::read_xlsx(here("plots", "99_paper_figs", "10c_snResolution_Top_Markers", "snResolution_top50MarkerGenes.xlsx"))
marker_stats[[1]] <- NULL
marker_stats |> count(cellType.target)

official_markers_tb <- tibble(cellType.short = names(official_markers),
                              cellType.target = names(official_markers),
                              gene = official_markers) |>
  unnest(gene) |>
  mutate(cellType.target = case_when(cellType.target == "Astro" ~'Astrocyte',
                                     cellType.target == "Micro" ~'Microglia',
                                     TRUE ~ cellType.target))

## Hb an Thal sub-type genes are datadriven (top mean ratio genes, glia is a mixed bag)
marker_details <- marker_stats |>
  right_join(official_markers_tb) |>
  arrange(rank_ratio) |>
  select(cellType.short, cellType.target, gene, rank_ratio, rank_marker) |>
  mutate(final_cell_type = cellType.target %in% marker_stats$cellType.target,
         anno = case_when(rank_ratio == 1 ~ "Data-Driven",
                          !final_cell_type | cellType.target == "Microglia" ~ "Literature",
                          TRUE ~ "Literature + Data-Supported"))

marker_details |> print(n = 22)

marker_stats |>
  group_by(cellType.target) |>
  arrange(rank_ratio) |>
  slice(1:5)

#### Mean Ratio Markers ####
load(here("processed-data", "06_deconvolution", "run_Bisque", "marker_stats_top_25_genes.Rdata"), verbose = TRUE)

marker_stats |> count(cellType.target)

markers_MR_top <- marker_stats |>
    filter(rank_ratio <= 5) |>
    select(cellType.target, gene, Symbol, rank_ratio) |>
    mutate(cellType.target = factor(cellType.target, levels = c("Astrocyte", "Endo","Microglia","Oligo","OPC",
                                                                "Excit.Thal", "Inhib.Thal", "LHb","MHb")),
           in_bulk = gene %in% rowData(rse_gene)$ensemblID,
           ) |>
    arrange(cellType.target, rank_ratio)

## some genes are missing
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
# [1] 39 68

marker_z_score <- scale(t(bulk_logcounts))

#### prep complex heatmap annotations  ####
# marker_method_colors <- c(`Data-Driven` = "#FFDA85",
#                           `Literature + Data-Supported` = "#F7A5A1",
#                           `Literature` = "#95E4EE")


## row sample annotations
sample_anno <- as.data.frame(colData(rse_gene), row.names = NULL) |>
    mutate(sn = BrNum %in% c("Br1092", "Br1204", "Br1469", "Br1735", "Br5555", "Br5558", "Br5639")) |>
    select(tot.Hb, PrimaryDx, sn)
rownames(sample_anno) <- rse_gene$BrNum

row_ha <- rowAnnotation(df = sample_anno)

## column marker gene annotations

broad_cell_type_colors <- c(sn_colors,   "MHb" = c("#F33A6A"), "LHb" = c("#0000FF"))[levels(markers_MR_top$cellType.target)]

column_ha <- HeatmapAnnotation(
    Marker_Gene = markers_MR_top$cellType.target,
    col = list(Marker_Gene = broad_cell_type_colors)
)

heatmapped <- Heatmap(marker_z_score,
                      name = "Z Score",
                      col = circlize::colorRamp2(seq(-4, 4,  8/10),
                                                 rev(RColorBrewer::brewer.pal(11, "RdBu"))),
                      # cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      right_annotation = row_ha,
                      top_annotation = column_ha,
                      column_split = markers_MR_top$cellType.target,
                      # heatmap_legend_param = list(
                      #   title = c("Z_Score"),
                      #   border = "black"
                      # )
                      )

# printing
pdf(here(plot_dir, "Hb_bulk_heatmap_MRtop5.pdf"), width = 12, height = 10)
  heatmapped
dev.off()


sessioninfo::session_info()
