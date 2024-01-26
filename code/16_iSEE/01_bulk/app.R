library("SummarizedExperiment")
library("iSEE")
library("shiny")

# ln -s ../../../processed-data/rse_objects/rse_gene_Habenula_Pilot.rda rse_gene_Habenula_Pilot.rda
load("rse_gene_Habenula_Pilot.rda", verbose = TRUE)

## Fix some phenotype data
rse_gene$PrimaryDx[rse_gene$PrimaryDx == "Schizo"] <- "SCZD"
rse_gene$Race[rse_gene$Race == "CAUC"] <- "EUA/CAUC"
rse_gene$AgeDeath <- round(rse_gene$AgeDeath, 2)

## Fix confusing names introduced at
## https://github.com/LieberInstitute/Habenula_Pilot/blob/d52c6e5f4fdf979d5d40589f54549691c5aac05e/code/02_bulk_qc/01_qc_data_Bukola.R#L41-L43
summary(rse_gene$numReads)
summary(rse_gene$numMapped)
colnames(colData(rse_gene))[colnames(colData(rse_gene)) == "numReads"] <- "log10NumReads"
colnames(colData(rse_gene))[colnames(colData(rse_gene)) == "numMapped"] <- "log10NumMapped"
rse_gene$numReads <- round(10^rse_gene$log10NumReads, 0)
rse_gene$numMapped <- round(10^rse_gene$log10NumMapped, 0)
summary(rse_gene$numReads)
summary(rse_gene$numMapped)

## iSEE configuration
initial <- list()

################################################################################
# Settings for Row data table 1
################################################################################

initial[["RowDataTable1"]] <- new("RowDataTable",
    Selected = "ENSG00000169676.5", Search = "",
    SearchColumns = c(
        "", "", "", "", "", "", "", "", "", "",
        ""
    ), HiddenColumns = c(
        "ensemblID", "gencodeID", "Class",
        "gencodeTx"
    ), VersionInfo = list(iSEE = structure(list(c(
        2L,
        12L, 0L
    )), class = c("package_version", "numeric_version"))), PanelId = c(RowDataTable = 1L), PanelHeight = 500L,
    PanelWidth = 6L, SelectionBoxOpen = FALSE, RowSelectionSource = "---",
    ColumnSelectionSource = "---", DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE,
    ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE,
    ColumnSelectionRestrict = FALSE, SelectionHistory = list()
)

################################################################################
# Settings for Feature assay plot 1
################################################################################

initial[["FeatureAssayPlot1"]] <- new("FeatureAssayPlot",
    Assay = "logcounts", XAxis = "Column data",
    XAxisColumnData = "PrimaryDx", XAxisFeatureName = "ENSG00000227232.5",
    XAxisFeatureSource = "---", XAxisFeatureDynamicSource = FALSE,
    YAxisFeatureName = "ENSG00000169676.5", YAxisFeatureSource = "RowDataTable1",
    YAxisFeatureDynamicSource = TRUE, FacetRowByColData = "BrNum",
    FacetColumnByColData = "BrNum", ColorByColumnData = "tot.Hb",
    ColorByFeatureNameAssay = "logcounts", ColorBySampleNameColor = "#FF0000",
    ShapeByColumnData = "BrNum", SizeByColumnData = "RIN", TooltipColumnData = character(0),
    FacetRowBy = "None", FacetColumnBy = "None", ColorBy = "Column data",
    ColorByDefaultColor = "#000000", ColorByFeatureName = "ENSG00000227232.5",
    ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE,
    ColorBySampleName = "R18346", ColorBySampleSource = "---",
    ColorBySampleDynamicSource = FALSE, ShapeBy = "None", SizeBy = "None",
    SelectionAlpha = 0.1, ZoomData = numeric(0), BrushData = list(),
    VisualBoxOpen = FALSE, VisualChoices = "Color", ContourAdd = FALSE,
    ContourColor = "#0000FF", PointSize = 1, PointAlpha = 1,
    Downsample = FALSE, DownsampleResolution = 200, CustomLabels = FALSE,
    CustomLabelsText = "R18346", FontSize = 1, LegendPointSize = 1,
    LegendPosition = "Bottom", HoverInfo = TRUE, LabelCenters = FALSE,
    LabelCentersBy = "BrNum", LabelCentersColor = "#000000",
    VersionInfo = list(iSEE = structure(list(c(2L, 12L, 0L)), class = c(
        "package_version",
        "numeric_version"
    ))), PanelId = c(FeatureAssayPlot = 1L),
    PanelHeight = 500L, PanelWidth = 6L, SelectionBoxOpen = FALSE,
    RowSelectionSource = "---", ColumnSelectionSource = "---",
    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE,
    RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE,
    SelectionHistory = list()
)

################################################################################
# Settings for Column data plot 1
################################################################################

initial[["ColumnDataPlot1"]] <- new("ColumnDataPlot",
    XAxis = "None", YAxis = "tot.Thal", XAxisColumnData = "BrNum",
    FacetRowByColData = "BrNum", FacetColumnByColData = "BrNum",
    ColorByColumnData = "Excit.Thal", ColorByFeatureNameAssay = "logcounts",
    ColorBySampleNameColor = "#FF0000", ShapeByColumnData = "BrNum",
    SizeByColumnData = "RIN", TooltipColumnData = character(0),
    FacetRowBy = "None", FacetColumnBy = "None", ColorBy = "Column data",
    ColorByDefaultColor = "#000000", ColorByFeatureName = "ENSG00000227232.5",
    ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE,
    ColorBySampleName = "R18346", ColorBySampleSource = "---",
    ColorBySampleDynamicSource = FALSE, ShapeBy = "None", SizeBy = "None",
    SelectionAlpha = 0.1, ZoomData = numeric(0), BrushData = list(),
    VisualBoxOpen = FALSE, VisualChoices = "Color", ContourAdd = FALSE,
    ContourColor = "#0000FF", PointSize = 1, PointAlpha = 1,
    Downsample = FALSE, DownsampleResolution = 200, CustomLabels = FALSE,
    CustomLabelsText = "R18346", FontSize = 1, LegendPointSize = 1,
    LegendPosition = "Bottom", HoverInfo = TRUE, LabelCenters = FALSE,
    LabelCentersBy = "BrNum", LabelCentersColor = "#000000",
    VersionInfo = list(iSEE = structure(list(c(2L, 12L, 0L)), class = c(
        "package_version",
        "numeric_version"
    ))), PanelId = c(ColumnDataPlot = 1L), PanelHeight = 500L,
    PanelWidth = 6L, SelectionBoxOpen = FALSE, RowSelectionSource = "---",
    ColumnSelectionSource = "ColumnDataTable1", DataBoxOpen = FALSE,
    RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE,
    RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE,
    SelectionHistory = list()
)

################################################################################
# Settings for Column data table 1
################################################################################

initial[["ColumnDataTable1"]] <- new("ColumnDataTable",
    Selected = "R18346", Search = "", SearchColumns = c(
        "",
        "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "",
        "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "",
        "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "",
        "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "",
        "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "",
        "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "",
        "", "", "", "", "", "", "", "", "", "", "", "", "", ""
    ), HiddenColumns = character(0),
    VersionInfo = list(iSEE = structure(list(c(2L, 12L, 0L)), class = c(
        "package_version",
        "numeric_version"
    ))), PanelId = c(ColumnDataTable = 1L),
    PanelHeight = 500L, PanelWidth = 6L, SelectionBoxOpen = FALSE,
    RowSelectionSource = "---", ColumnSelectionSource = "---",
    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE,
    RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE,
    SelectionHistory = list()
)

## Build the iSEE app
iSEE(rse_gene, appTitle = "habenulaPilot - bulk RNA-seq", initial = initial)
