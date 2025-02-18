library("SingleCellExperiment")
library("iSEE")
library("shiny")

## Load data
sce <- readRDS("sce.rds")
sn_colors <- readRDS("sn_colors.rds")
bulk_colors <- readRDS("bulk_colors.rds")

## iSEE configuration
initial <- list()

################################################################################
# Settings for Reduced dimension plot 1
################################################################################

initial[["ReducedDimensionPlot1"]] <- new("ReducedDimensionPlot", Type = "PCA", XAxis = 2L, YAxis = 3L,
    FacetRowByColData = "Sample", FacetColumnByColData = "Sample",
    ColorByColumnData = "final_Annotations", ColorByFeatureNameAssay = "logcounts",
    ColorBySampleNameColor = "#FF0000", ShapeByColumnData = "Sample",
    SizeByColumnData = "sum", TooltipColumnData = character(0),
    FacetRowBy = "None", FacetColumnBy = "None", ColorBy = "Column data",
    ColorByDefaultColor = "#000000", ColorByFeatureName = "MIR1302-2HG",
    ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE,
    ColorBySampleName = "Br1092_AAACCCAAGTCTACCA-1", ColorBySampleSource = "---",
    ColorBySampleDynamicSource = FALSE, ShapeBy = "None", SizeBy = "None",
    SelectionAlpha = 0.1, ZoomData = numeric(0), BrushData = list(
        xmin = 97.971656170844, xmax = 112.24228661805, ymin = -28.754901934858,
        ymax = 37.770993498102, coords_css = list(xmin = 683.199999809265,
            xmax = 720.199999809265, ymin = 129.799995422363,
            ymax = 215.799995422363), coords_img = list(xmin = 1707.99999952316,
            xmax = 1800.49999952316, ymin = 324.499988555908,
            ymax = 539.499988555908), img_css_ratio = list(x = 2.5,
            y = 2.5), mapping = list(x = "X", y = "Y", colour = "ColorBy"),
        domain = list(left = -148.420638687221, right = 214.343073489101,
            bottom = -139.50748746747, top = 119.180968334431),
        range = list(left = 110.923587328767, right = 2462.30136986301,
            bottom = 897.432876712329, top = 61.3972602739727),
        log = list(x = NULL, y = NULL), direction = "xy", brushId = "ReducedDimensionPlot1_Brush",
        outputId = "ReducedDimensionPlot1"), VisualBoxOpen = FALSE,
    VisualChoices = "Color", ContourAdd = FALSE, ContourColor = "#0000FF",
    FixAspectRatio = FALSE, ViolinAdd = TRUE, PointSize = 1,
    PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200,
    CustomLabels = FALSE, CustomLabelsText = "Br1092_AAACCCAAGTCTACCA-1",
    FontSize = 1, LegendPointSize = 1, LegendPosition = "Bottom",
    HoverInfo = TRUE, LabelCenters = FALSE, LabelCentersBy = "Sample",
    LabelCentersColor = "#000000", VersionInfo = list(iSEE = structure(list(
        c(2L, 18L, 0L)), class = c("package_version", "numeric_version"
    ))), PanelId = c(ReducedDimensionPlot = 1L), PanelHeight = 500L,
    PanelWidth = 6L, SelectionBoxOpen = FALSE, RowSelectionSource = "---",
    ColumnSelectionSource = "---", DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE,
    ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE,
    ColumnSelectionRestrict = FALSE, SelectionHistory = list())

################################################################################
# Settings for Complex heatmap 1
################################################################################

initial[["ComplexHeatmapPlot1"]] <- new("ComplexHeatmapPlot", Assay = "logcounts", CustomRows = TRUE,
    CustomRowsText = "POU4F1\nGPR151\nCHRNB4\nHTR2C\nLYPD6B\nADARB2\nRORB\nSYT1\nSLC17A6\nGAD1\nMOBP\nPDGFRA\nAQP4\nITIH5\nCSF1R",
    ClusterRows = TRUE, ClusterRowsDistance = "spearman", ClusterRowsMethod = "ward.D2",
    DataBoxOpen = FALSE, VisualChoices = "Annotations", ColumnData = "final_Annotations",
    RowData = character(0), CustomBounds = FALSE, LowerBound = NA_real_,
    UpperBound = NA_real_, AssayCenterRows = TRUE, AssayScaleRows = TRUE,
    DivergentColormap = "purple < black < yellow", ShowDimNames = "Rows",
    LegendPosition = "Right", LegendDirection = "Vertical", VisualBoxOpen = FALSE,
    NamesRowFontSize = 10, NamesColumnFontSize = 10, ShowColumnSelection = TRUE,
    OrderColumnSelection = TRUE, VersionInfo = list(iSEE = structure(list(
        c(2L, 18L, 0L)), class = c("package_version", "numeric_version"
    ))), PanelId = c(ComplexHeatmapPlot = 1L), PanelHeight = 500L,
    PanelWidth = 6L, SelectionBoxOpen = FALSE, RowSelectionSource = "---",
    ColumnSelectionSource = "---", RowSelectionDynamicSource = FALSE,
    ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE,
    ColumnSelectionRestrict = FALSE, SelectionHistory = list())

################################################################################
# Settings for Row data table 1
################################################################################

initial[["RowDataTable1"]] <- new("RowDataTable", Selected = "CHAT", Search = "", SearchColumns = c("",
"", "", "", "", "", "", ""), HiddenColumns = "Type", VersionInfo = list(
    iSEE = structure(list(c(2L, 18L, 0L)), class = c("package_version",
    "numeric_version"))), PanelId = c(RowDataTable = 1L), PanelHeight = 500L,
    PanelWidth = 6L, SelectionBoxOpen = FALSE, RowSelectionSource = "---",
    ColumnSelectionSource = "---", DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE,
    ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE,
    ColumnSelectionRestrict = FALSE, SelectionHistory = list())

################################################################################
# Settings for Feature assay plot 1
################################################################################

initial[["FeatureAssayPlot1"]] <- new("FeatureAssayPlot", Assay = "logcounts", XAxis = "Column data",
    XAxisColumnData = "final_Annotations", XAxisFeatureName = "MIR1302-2HG",
    XAxisFeatureSource = "---", XAxisFeatureDynamicSource = FALSE,
    YAxisFeatureName = "CHAT", YAxisFeatureSource = "RowDataTable1",
    YAxisFeatureDynamicSource = FALSE, FacetRowByColData = "Sample",
    FacetColumnByColData = "Sample", ColorByColumnData = "final_Annotations",
    ColorByFeatureNameAssay = "logcounts", ColorBySampleNameColor = "#FF0000",
    ShapeByColumnData = "Sample", SizeByColumnData = "sum", TooltipColumnData = character(0),
    FacetRowBy = "None", FacetColumnBy = "None", ColorBy = "Column data",
    ColorByDefaultColor = "#000000", ColorByFeatureName = "MIR1302-2HG",
    ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE,
    ColorBySampleName = "Br1092_AAACCCAAGTCTACCA-1", ColorBySampleSource = "---",
    ColorBySampleDynamicSource = FALSE, ShapeBy = "None", SizeBy = "None",
    SelectionAlpha = 0.1, ZoomData = numeric(0), BrushData = list(),
    VisualBoxOpen = FALSE, VisualChoices = "Color", ContourAdd = FALSE,
    ContourColor = "#0000FF", FixAspectRatio = FALSE, ViolinAdd = TRUE,
    PointSize = 1, PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200,
    CustomLabels = FALSE, CustomLabelsText = "Br1092_AAACCCAAGTCTACCA-1",
    FontSize = 1, LegendPointSize = 1, LegendPosition = "Bottom",
    HoverInfo = TRUE, LabelCenters = FALSE, LabelCentersBy = "Sample",
    LabelCentersColor = "#000000", VersionInfo = list(iSEE = structure(list(
        c(2L, 18L, 0L)), class = c("package_version", "numeric_version"
    ))), PanelId = c(FeatureAssayPlot = 1L), PanelHeight = 500L,
    PanelWidth = 6L, SelectionBoxOpen = FALSE, RowSelectionSource = "---",
    ColumnSelectionSource = "---", DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE,
    ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE,
    ColumnSelectionRestrict = FALSE, SelectionHistory = list())

## Build the iSEE app
iSEE(sce,
        appTitle = "habenulaPilot - snRNA-seq", initial = initial,
        colormap = ExperimentColorMap(
            colData = list(
                final_Annotations = function(n) {
                    return(sn_colors)
                },
                broad_Annotations = function(n) {
                    return(bulk_colors)
                }
            )
        )
    )
