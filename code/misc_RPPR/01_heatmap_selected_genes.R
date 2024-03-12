library("here")
library("SingleCellExperiment")
library("tidyverse")
library("ComplexHeatmap")
library("sessioninfo")

# creating plot directory
plot_dir <- here("plots", "misc_RPPR")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)


# sourcing official color palette
source(file = here("code", "99_paper_figs", "source_colors.R"))

# Load pseudobulked data
load(
    here(
        "processed-data",
        "04_snRNA-seq",
        "sce_objects",
        "sce_pseudobulk_final_Annotations.Rdata"
    ),
    verbose = TRUE
)

# list of marker genes
# official_markers <- list(
#     "Oligo" = c("MOBP"),
#     "OPC" = c("PDGFRA"),
#     "Micro" = c("CSF1R"),
#     "Astro" = c("AQP4"),
#     "Endo" = c("ITIH5"),
#     "Thal" = c("LYPD6B", "ADARB2", "RORB"),
#     # "LHb.1" = c("LINC02653"), #  , ATP8B1
#     # "LHb.2" = c("AC073071.1"),
#     # "LHb.3" = c ("ENTHD1"),
#     # "LHb.4" = c("TLE2"),
#     # "LHb.5" = c("LINC01619"),
#     # "LHb.6" = c("TACR3"),
#     # "LHb.7" = c("AC008619.1"),
#     # "MHb.1" = c("EXOC1L"),
#     # "MHb.2" = c("CHAT"),
#     # "MHb.3" = c("BHLHE22"),
#     "Hb" = c("POU4F1", "GPR151", "TAC3"),
#     # BARHL1
#     "MHb" = c("CHRNB4"),
#     "LHb" = c("HTR2C"),
#     "Neuron" = c("SYT1"),
#     "Exc_Neuron" = c("SLC17A6"),
#     "Inh_Neuron" = c("GAD1")
# )

official_markers  <- list(
    "Hb" = c("MMRN1", "GPR151", "POU4F1"),
    "LHb" = c("EPHA5", "TLL1", "AC109466.1", "AC008415.1", "GPR149", "GNG8", "NEUROD1"),
    "MHb" = c("RASGRP1", "SLC5A7", "CHRNB3", "SCUBE1", "LINC02143", "CD24")
    # ,
    # "Oligo" = c("MOBP"),
    # "OPC" = c("PDGFRA"),
    # "Micro" = c("CSF1R"),
    # "Astro" = c("AQP4"),
    # "Endo" = c("ITIH5"),
    # "Thal" = c("LYPD6B", "ADARB2", "RORB"),
    # "Neuron" = c("SYT1"),
    # "Exc_Neuron" = c("SLC17A6"),
    # "Inh_Neuron" = c("GAD1")
)

#### Check marker genes ####
marker_stats <-
    readxl::read_xlsx(
        here(
            "plots",
            "99_paper_figs",
            "10c_snResolution_Top_Markers",
            "snResolution_top50MarkerGenes.xlsx"
        )
    )
marker_stats[[1]] <- NULL
marker_stats |> count(cellType.target)

official_markers_tb <-
    tibble(
        cellType.short = names(official_markers),
        cellType.target = names(official_markers),
        gene = official_markers
    ) |>
    unnest(gene) |>
    mutate(
        cellType.target = case_when(
            cellType.target == "Astro" ~ "Astrocyte",
            cellType.target == "Micro" ~ "Microglia",
            TRUE ~ cellType.target
        )
    )

## Hb an Thal sub-type genes are datadriven (top mean ratio genes, glia is a mixed bag)
marker_details <- marker_stats |>
    right_join(official_markers_tb) |>
    arrange(rank_ratio) |>
    select(cellType.short, cellType.target, gene, rank_ratio, rank_marker) |>
    mutate(
        final_cell_type = cellType.target %in% marker_stats$cellType.target,
        anno = case_when(
            rank_ratio == 1 ~ "Data-Driven",
            !final_cell_type |
                cellType.target == "Microglia" ~ "Literature",
            TRUE ~ "Literature + Data-Supported"
        )
    )

marker_details |> print(n = 22)

marker_stats |>
    group_by(cellType.target) |>
    arrange(rank_ratio) |>
    slice(1:5)

#### prep complex heatmap annotations  ####
# explicit color scheme
color_official_markers <- c(
    "Oligo" = c("#4d5802"),
    "OPC" = c("#d3c871"),
    "Micro" = c("#222222"),
    "Astro" = c("#8d363c"),
    "Endo" = c("#ee6c14"),
    "Thal" = c("#EADDCA"),
    "LHb.1" = c("#0085af"),
    "LHb.2" = c("#0096FF"),
    "LHb.3" = c("#89CFF0"),
    "LHb.4" = c("#6F8FAF"),
    "LHb.5" = c("#40E0D0"),
    "LHb.6" = c("#008080"),
    "LHb.7" = c("#7DF9FF"),
    "MHb.1" = c("#FF00FF"),
    "MHb.2" = c("#FAA0A0"),
    "MHb.3" = c("#fa246a"),
    "Hb" = c("#702963"),
    "MHb" = c("#F33A6A"),
    "LHb" = c("#0000FF"),
    "Neuron" = c("#5C5C5C"),
    "Exc_Neuron" = c("#8F8F8F"),
    "Inh_Neuron" = c("#C2C2C2")
)

marker_method_colors <- c(
    `Data-Driven` = "#FFDA85",
    `Literature + Data-Supported` = "#F7A5A1",
    `Literature` = "#95E4EE"
)

####### PLOTTING ###############################################################
# Plotting ComplexHeatmap
sce <- sce_pb
clusterMethod <- "final_Annotations"
markerList <- official_markers

# Replacing genes with symbols for heatmap (remember, this is pseudobulked data)
rownames(sce) <- rowData(sce)$Symbol

# renaming rownnames of colData(sce) based on row annotations
rownames(colData(sce)) <- colData(sce)[, clusterMethod]

# Making data frame of genes we are interested in annd their general classification
# markTable <- as.data.frame(unlist(markerList)) |>
#   rownames_to_column("cellType") |>
#   rename(gene = `unlist(markerList)`) |>
#   # mutate(cellType = gsub("\\d+", "", cellType)) |>
#   filter(gene %in% rowData(sce_reorder)$Symbol)

markTable <- marker_details |>
    mutate(
        cellType.short = factor(cellType.short, levels = names(official_markers)),
        anno = factor(
            anno,
            levels = c("Data-Driven", "Literature + Data-Supported", "Literature")
        )
    ) |>
    arrange(anno, cellType.short)
## fix order in colsplit
official_markers_order <- c(
    "Hb",
    "LHb",
    "MHb",
    "Thal",
    "Neuron",
    "Exc_Neuron",
    "Inh_Neuron",
    "Oligo",
    "OPC",
    "Astro",
    "Endo",
    "Micro"
)
# stopifnot(all(official_markers_order %in% names(official_markers)))
markTable$cellType.short <- factor(markTable$cellType.short,
    levels = official_markers_order
)

markTable |> print(n = 22)

row_namers <- c(
    "LHb.1",
    "LHb.2",
    "LHb.3",
    "LHb.4",
    "LHb.5",
    "LHb.6",
    "LHb.7",
    "MHb.1",
    "MHb.2",
    "MHb.3",
    "Inhib.Thal",
    "Excit.Thal",
    "Oligo",
    "OPC",
    "Astrocyte",
    "Endo",
    "Microglia"
)

# row_namers <- markTable |> filter(final_cell_type) |> pull(cellType.target)
# marker_order <- markTable |> pull(gene)

# reordering sce object for plottability
sce_reorder <- sce[markTable$gene, row_namers]
sce_reorder$final_Annotations

# getting z scores
marker_z_score <- scale(t(logcounts(sce_reorder)))
# corner(marker_z_score)

identical(markTable$gene, colnames(marker_z_score))

# heatmap columns annotation
column_ha <- HeatmapAnnotation(
    Marker_Gene = markTable$cellType.short,
    # Marker_Method = markTable$anno,
    col = list(
        Marker_Gene = color_official_markers,
        Marker_Method = marker_method_colors
    )
)

# grabbing the annotations per cluster from the sce_reorder object
# clusterData <- as.data.frame(colData(sce_reorder)[,clusterMethod])
# names(clusterData) <- "cellType"

# prepping the colors we want
# for cell type
# num_pal <- length(unique(clusterData$cellType))
# col_pal_ct <- grabColors(num_pal, start = 4)
# names(col_pal_ct) = unique(clusterData$cellType)

# heatmap row annotationn
# identical(clusterData$cellType, rownames(marker_z_score))

row_ha <- rowAnnotation(
    Clusters = rownames(marker_z_score),
    col = list(Clusters = sn_colors)
)

heatmapped <- Heatmap(
    marker_z_score,
    name = "Z Score",
    col = circlize::colorRamp2(
        seq(-4, 4, 8 / 10),
        rev(RColorBrewer::brewer.pal(11, "RdBu"))
    ),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    right_annotation = row_ha,
    top_annotation = column_ha,
    column_split = markTable$cellType.short,
    column_title_rot = 30
    # heatmap_legend_param = list(
    #   title = c("Z_Score"),
    #   border = "black"
    # )
)

# printing
pdf(
    here(plot_dir, "RPPR_2024_heatmap.pdf"),
    width = 12,
    height = 8
)
heatmapped
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.3.2 (2023-10-31)
#  os       macOS Sonoma 14.3.1
#  system   aarch64, darwin20
#  ui       RStudio
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       America/New_York
#  date     2024-03-12
#  rstudio  2023.12.1+402 Ocean Storm (desktop)
#  pandoc   3.1.12.1 @ /opt/homebrew/bin/pandoc
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version    date (UTC) lib source
#  abind                  1.4-5      2016-07-21 [1] CRAN (R 4.3.0)
#  Biobase              * 2.62.0     2023-10-26 [1] Bioconductor
#  BiocGenerics         * 0.48.1     2023-11-02 [1] Bioconductor
#  biocthis               1.12.0     2023-10-26 [1] Bioconductor
#  bitops                 1.0-7      2021-04-24 [1] CRAN (R 4.3.0)
#  brio                   1.1.4      2023-12-10 [1] CRAN (R 4.3.1)
#  cachem                 1.0.8      2023-05-01 [1] CRAN (R 4.3.0)
#  Cairo                  1.6-2      2023-11-28 [1] CRAN (R 4.3.1)
#  cellranger             1.1.0      2016-07-27 [1] CRAN (R 4.3.0)
#  circlize               0.4.16     2024-02-20 [1] CRAN (R 4.3.1)
#  cli                    3.6.2      2023-12-11 [1] CRAN (R 4.3.1)
#  clue                   0.3-65     2023-09-23 [1] CRAN (R 4.3.1)
#  cluster                2.1.6      2023-12-01 [1] CRAN (R 4.3.1)
#  codetools              0.2-19     2023-02-01 [1] CRAN (R 4.3.2)
#  colorout             * 1.3-0.2    2024-02-27 [1] Github (jalvesaq/colorout@c6113a2)
#  colorspace             2.1-0      2023-01-23 [1] CRAN (R 4.3.0)
#  ComplexHeatmap       * 2.18.0     2023-10-26 [1] Bioconductor
#  crayon                 1.5.2      2022-09-29 [1] CRAN (R 4.3.0)
#  data.table             1.15.0     2024-01-30 [1] CRAN (R 4.3.1)
#  DelayedArray           0.28.0     2023-11-06 [1] Bioconductor
#  devtools             * 2.4.5      2022-10-11 [1] CRAN (R 4.3.0)
#  digest                 0.6.34     2024-01-11 [1] CRAN (R 4.3.1)
#  doParallel             1.0.17     2022-02-07 [1] CRAN (R 4.3.0)
#  dplyr                * 1.1.4      2023-11-17 [1] CRAN (R 4.3.1)
#  ellipsis               0.3.2      2021-04-29 [1] CRAN (R 4.3.0)
#  fansi                  1.0.6      2023-12-08 [1] CRAN (R 4.3.1)
#  fastmap                1.1.1      2023-02-24 [1] CRAN (R 4.3.0)
#  forcats              * 1.0.0      2023-01-29 [1] CRAN (R 4.3.0)
#  foreach                1.5.2      2022-02-02 [1] CRAN (R 4.3.0)
#  fs                     1.6.3      2023-07-20 [1] CRAN (R 4.3.0)
#  generics               0.1.3      2022-07-05 [1] CRAN (R 4.3.0)
#  GenomeInfoDb         * 1.38.6     2024-02-10 [1] Bioconductor 3.18 (R 4.3.2)
#  GenomeInfoDbData       1.2.11     2024-02-27 [1] Bioconductor
#  GenomicRanges        * 1.54.1     2023-10-30 [1] Bioconductor
#  GetoptLong             1.0.5      2020-12-15 [1] CRAN (R 4.3.0)
#  ggplot2              * 3.5.0      2024-02-23 [1] CRAN (R 4.3.1)
#  GlobalOptions          0.1.2      2020-06-10 [1] CRAN (R 4.3.0)
#  glue                   1.7.0      2024-01-09 [1] CRAN (R 4.3.1)
#  gtable                 0.3.4      2023-08-21 [1] CRAN (R 4.3.0)
#  here                 * 1.0.1      2020-12-13 [1] CRAN (R 4.3.0)
#  hms                    1.1.3      2023-03-21 [1] CRAN (R 4.3.0)
#  htmltools              0.5.7      2023-11-03 [1] CRAN (R 4.3.1)
#  htmlwidgets            1.6.4      2023-12-06 [1] CRAN (R 4.3.1)
#  httpuv                 1.6.14     2024-01-26 [1] CRAN (R 4.3.1)
#  IRanges              * 2.36.0     2023-10-26 [1] Bioconductor
#  iterators              1.0.14     2022-02-05 [1] CRAN (R 4.3.0)
#  later                  1.3.2      2023-12-06 [1] CRAN (R 4.3.1)
#  lattice                0.22-5     2023-10-24 [1] CRAN (R 4.3.1)
#  lifecycle              1.0.4      2023-11-07 [1] CRAN (R 4.3.1)
#  lubridate            * 1.9.3      2023-09-27 [1] CRAN (R 4.3.1)
#  magick                 2.8.3      2024-02-18 [1] CRAN (R 4.3.1)
#  magrittr               2.0.3      2022-03-30 [1] CRAN (R 4.3.0)
#  Matrix                 1.6-5      2024-01-11 [1] CRAN (R 4.3.2)
#  MatrixGenerics       * 1.14.0     2023-10-26 [1] Bioconductor
#  matrixStats          * 1.2.0      2023-12-11 [1] CRAN (R 4.3.1)
#  memoise                2.0.1      2021-11-26 [1] CRAN (R 4.3.0)
#  mime                   0.12       2021-09-28 [1] CRAN (R 4.3.0)
#  miniUI                 0.1.1.1    2018-05-18 [1] CRAN (R 4.3.0)
#  munsell                0.5.0      2018-06-12 [1] CRAN (R 4.3.0)
#  pillar                 1.9.0      2023-03-22 [1] CRAN (R 4.3.0)
#  pkgbuild               1.4.3      2023-12-10 [1] CRAN (R 4.3.1)
#  pkgconfig              2.0.3      2019-09-22 [1] CRAN (R 4.3.0)
#  pkgload                1.3.4      2024-01-16 [1] CRAN (R 4.3.1)
#  png                    0.1-8      2022-11-29 [1] CRAN (R 4.3.0)
#  profvis                0.3.8      2023-05-02 [1] CRAN (R 4.3.0)
#  promises               1.2.1      2023-08-10 [1] CRAN (R 4.3.0)
#  prompt                 1.0.2.9000 2024-02-27 [1] Github (gaborcsardi/prompt@17bd0e1)
#  purrr                * 1.0.2      2023-08-10 [1] CRAN (R 4.3.0)
#  R.cache                0.16.0     2022-07-21 [1] CRAN (R 4.3.0)
#  R.methodsS3            1.8.2      2022-06-13 [1] CRAN (R 4.3.0)
#  R.oo                   1.26.0     2024-01-24 [1] CRAN (R 4.3.1)
#  R.utils                2.12.3     2023-11-18 [1] CRAN (R 4.3.1)
#  R6                     2.5.1      2021-08-19 [1] CRAN (R 4.3.0)
#  RColorBrewer           1.1-3      2022-04-03 [1] CRAN (R 4.3.0)
#  Rcpp                   1.0.12     2024-01-09 [1] CRAN (R 4.3.1)
#  RCurl                  1.98-1.14  2024-01-09 [1] CRAN (R 4.3.1)
#  readr                * 2.1.5      2024-01-10 [1] CRAN (R 4.3.1)
#  readxl                 1.4.3      2023-07-06 [1] CRAN (R 4.3.0)
#  remotes                2.4.2.1    2023-07-18 [1] CRAN (R 4.3.0)
#  rjson                  0.2.21     2022-01-09 [1] CRAN (R 4.3.0)
#  rlang                  1.1.3      2024-01-10 [1] CRAN (R 4.3.1)
#  rprojroot              2.0.4      2023-11-05 [1] CRAN (R 4.3.1)
#  rsthemes               0.4.0      2024-02-27 [1] Github (gadenbuie/rsthemes@34a55a4)
#  rstudioapi             0.15.0     2023-07-07 [1] CRAN (R 4.3.0)
#  S4Arrays               1.2.0      2023-10-26 [1] Bioconductor
#  S4Vectors            * 0.40.2     2023-11-25 [1] Bioconductor 3.18 (R 4.3.2)
#  scales                 1.3.0      2023-11-28 [1] CRAN (R 4.3.1)
#  sessioninfo          * 1.2.2      2021-12-06 [1] CRAN (R 4.3.0)
#  shape                  1.4.6.1    2024-02-23 [1] CRAN (R 4.3.1)
#  shiny                  1.8.0      2023-11-17 [1] CRAN (R 4.3.1)
#  SingleCellExperiment * 1.24.0     2023-11-06 [1] Bioconductor
#  SparseArray            1.2.4      2024-02-10 [1] Bioconductor 3.18 (R 4.3.2)
#  stringi                1.8.3      2023-12-11 [1] CRAN (R 4.3.1)
#  stringr              * 1.5.1      2023-11-14 [1] CRAN (R 4.3.1)
#  styler                 1.10.2     2023-08-29 [1] CRAN (R 4.3.0)
#  SummarizedExperiment * 1.32.0     2023-11-06 [1] Bioconductor
#  suncalc                0.5.1      2022-09-29 [1] CRAN (R 4.3.0)
#  testthat             * 3.2.1      2023-12-02 [1] CRAN (R 4.3.1)
#  tibble               * 3.2.1      2023-03-20 [1] CRAN (R 4.3.0)
#  tidyr                * 1.3.1      2024-01-24 [1] CRAN (R 4.3.1)
#  tidyselect             1.2.0      2022-10-10 [1] CRAN (R 4.3.0)
#  tidyverse            * 2.0.0      2023-02-22 [1] CRAN (R 4.3.0)
#  timechange             0.3.0      2024-01-18 [1] CRAN (R 4.3.1)
#  tzdb                   0.4.0      2023-05-12 [1] CRAN (R 4.3.0)
#  urlchecker             1.0.1      2021-11-30 [1] CRAN (R 4.3.0)
#  usethis              * 2.2.3      2024-02-19 [1] CRAN (R 4.3.1)
#  utf8                   1.2.4      2023-10-22 [1] CRAN (R 4.3.1)
#  vctrs                  0.6.5      2023-12-01 [1] CRAN (R 4.3.1)
#  withr                  3.0.0      2024-01-16 [1] CRAN (R 4.3.1)
#  xtable                 1.8-4      2019-04-21 [1] CRAN (R 4.3.0)
#  XVector                0.42.0     2023-10-26 [1] Bioconductor
#  zlibbioc               1.48.0     2023-10-26 [1] Bioconductor
#
#  [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
