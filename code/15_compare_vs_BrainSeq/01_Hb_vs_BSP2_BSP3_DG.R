library("here")
library("BiocFileCache")
library("readxl")
library("dplyr")
library("GGally")
library("ggrepel")
library("patchwork")
library("sessioninfo")

## Output directories
dir_rdata <- here("processed-data", "15_compare_vs_BrainSeq")
dir_plots <- here("plots", "15_compare_vs_BrainSeq")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

## Inputs
dir_input <- here("processed-data", "10_DEA", "04_DEA")
# dir_rawdata <- "/dcs04/lieber/lcolladotor/qSVA_LIBD3080/qsva_brain/brainseq_phase2_qsv/rdas"
dir_rawdata <- here("raw-data", "15_compare_vs_BrainSeq")

## Read in Habenula data: DLPFC and HIPPO (Hippocampus)
habenula <- read.table(file.path(dir_input, "DEA_All-gene_qc-totAGene-qSVs-Hb-Thal.tsv"), header = TRUE)

## Read in BSP2 data
load(file.path(dir_rawdata, "dxStats_dlpfc_filtered_qSVA_geneLevel_noHGoldQSV_matchDLPFC.rda"), verbose = TRUE)
bsp2_dlpfc <- outGene
rm(outGene, outGene0, outGeneNoAdj)
load(file.path(dir_rawdata, "dxStats_hippo_filtered_qSVA_geneLevel_noHGoldQSV_matchHIPPO.rda"), verbose = TRUE)
bsp2_hpc <- outGene
rm(outGene, outGene0, outGeneNoAdj)

## Read in BSP3 data: Caudate
## Data URL from Table S10 at https://erwinpaquolalab.libd.org/caudate_eqtl/
bfc <- BiocFileCache()
bsp3_caudate_file <- BiocFileCache::bfcrpath(
    bfc,
    "https://caudate-eqtl.s3.us-west-2.amazonaws.com/BrainSeq_Phase3_Caudate_DifferentialExpression_DxSZ_all.txt.gz",
    exact = TRUE
)
bsp3_caudate <- read.table(bsp3_caudate_file, header = TRUE, fill = TRUE)
bsp3_caudate <- bsp3_caudate[bsp3_caudate$Type == "Gene", ]

## DG-GCL results
## Table S10 from https://www.nature.com/articles/s41593-020-0604-z#MOESM2
dg <- read_xlsx(file.path(dir_rawdata, "41593_2020_604_MOESM2_ESM.xlsx"), sheet = "Table_S10")


## Function to build cor table in long format
cor_long <- function(input) {
    x <- cor(input,
        use = "pairwise.complete.obs",
        method = "pearson")
    rowCol <- expand.grid(rownames(x), colnames(x))
    labs <- rowCol[as.vector(upper.tri(x, diag = FALSE)), ]
    df <- cbind(labs, x[upper.tri(x, diag = FALSE)])
    colnames(df) <- c("Row", "Col", "Cor")
    return(df)
}

## Compute it by subsets
cor_long_type <- function(all_input) {
    regions <- c("Hb", "DLPFC", "HIPPO", "Caudate", "DG")
    res <- rename(cor_long(all_input[, regions]), All = Cor)
    res <- left_join(
        res,
        rename(cor_long(all_input[all_input$adj.P.Val < 0.05, regions]), Sig = Cor),
        by = c("Row", "Col")
    )
    left_join(
        res,
        rename(cor_long(all_input[all_input$adj.P.Val >= 0.05, regions]), NotSig = Cor),
        by = c("Row", "Col")
    )
}

## Region colors from https://github.com/LieberInstitute/Habenula_Pilot/blob/dd21b1b6b0855d4688d625121c1e1de7a3836b18/code/03_bulk_pca/02_multiregion_PCA.R#L140C1-L152C21
region_colors <- c(Amygdala = "#ff9ccb",
    BLA = "#c10040",
    CA = "#ff7168",
    MeA = "#b9008b",
    DG = "#00960e",
    HIPPO = "#99C71A",
    dACC = "#0094fc",
    sACC = "#014abf",
    DLPFC = "#c495ff",
    mPFC = "#8330b6",
    Caudate = "#65717B",
    Hb = "#F4D23E" #Mustard
)

## Make a plot
plot_cor_type <- function(all_input, hb_only = TRUE) {
    df <- cor_long_type(all_input)
    if(hb_only) {
        df <- subset(df, Row == "Hb")
    }
    t_res <- t.test(df[, c("Sig", "NotSig")])
    # print(t_res)
    lims <- range(c(df$Sig, df$NotSig))
    ggplot(df, aes(x = NotSig, y = Sig, colour = Col)) +
        geom_point(size = 5) +
        xlab("cor with Hb FDR >= 0.05") +
        ylab("cor with Hb FDR < 0.05") +
        scale_color_manual(values = region_colors) +
        guides(colour = guide_legend(title="Region")) +
        geom_abline(intercept = 0, slope = 1, col = "red") +
        xlim(lims) +
        ylim(lims) +
        annotate("text", label = paste0("p = ", signif(t_res$p.value, 3)), x = max(lims) - 0.08 * max(lims), y = min(lims), size = 5) +
        theme_bw(base_size = 20)
}


## Merge data: FDRs
all_FDR <- mutate(habenula, Hb = adj.P.Val)[, c("gencodeID", "Symbol", "Hb")]
all_FDR <- left_join(all_FDR, mutate(bsp2_dlpfc, DLPFC = adj.P.Val)[, c("gencodeID", "DLPFC")], by = "gencodeID")
all_FDR <- left_join(all_FDR, mutate(bsp2_hpc, HIPPO = adj.P.Val)[, c("gencodeID", "HIPPO")], by = "gencodeID")
all_FDR <- left_join(all_FDR, mutate(bsp3_caudate, Caudate = as.numeric(adj.P.Val))[, c("gencodeID", "Caudate")], by = "gencodeID")
all_FDR <- left_join(all_FDR, mutate(dg, DG = SZ_adj.P.Val)[, c("gencodeID", "DG")], by = "gencodeID")

addmargins(table(all_FDR$Hb < 0.05, all_FDR$DLPFC < 0.05, useNA = "ifany"))
#       FALSE  TRUE  <NA>   Sum
# FALSE 20452   235  2024 22711
# TRUE     37     1     7    45
# Sum   20489   236  2031 22756
all_FDR[which(all_FDR$Hb < 0.05 & all_FDR$DLPFC < 0.05), ]
#                gencodeID Symbol         Hb      DLPFC      HIPPO   Caudate        DG
# 20138 ENSG00000117877.10 CD3EAP 0.04304279 0.04036771 0.03465433 0.5060051 0.4006797


addmargins(table(all_FDR$Hb < 0.05, all_FDR$HIPPO < 0.05, useNA = "ifany"))
#       FALSE  TRUE  <NA>   Sum
# FALSE 20642    45  2024 22711
# TRUE     37     1     7    45
# Sum   20679    46  2031 22756
all_FDR[which(all_FDR$Hb < 0.05 & all_FDR$HIPPO < 0.05), ]
#                gencodeID Symbol         Hb      DLPFC      HIPPO   Caudate        DG
# 20138 ENSG00000117877.10 CD3EAP 0.04304279 0.04036771 0.03465433 0.5060051 0.4006797

## More info on this gene:
# https://www.genecards.org/cgi-bin/carddisp.pl?gene=POLR1G&keywords=ENSG00000117877
# https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000117877;r=19:45406644-45410737
#
# ## There are no hits on PubMed for this gene with SCZD
# https://pubmed.ncbi.nlm.nih.gov/?term=CD3EAP+schizophrenia
# https://pubmed.ncbi.nlm.nih.gov/?term=POLR1G+schizophrenia


addmargins(table(all_FDR$Hb < 0.05, all_FDR$Caudate < 0.05, useNA = "ifany"))
#       FALSE  TRUE  <NA>   Sum
# FALSE 15297  2425  4989 22711
# TRUE     24    12     9    45
# Sum   15321  2437  4998 22756
all_FDR[which(all_FDR$Hb < 0.05 & all_FDR$Caudate < 0.05), ]
#                gencodeID   Symbol         Hb     DLPFC     HIPPO      Caudate        DG
# 313   ENSG00000169641.13    LUZP1 0.03138251 0.1679386 0.8085046 2.779692e-02 0.9904700
# 659    ENSG00000186603.5     HPDL 0.04304279        NA        NA 2.568891e-02        NA
# 3142  ENSG00000153250.18    RBMS1 0.04304279 0.8305017 0.7584133 6.451821e-04        NA
# 6697  ENSG00000070614.14    NDST1 0.04304279 0.9697855 0.9895913 4.947168e-02 0.9312463
# 10758  ENSG00000160447.6     PKN3 0.04304279 0.6026409 0.9441769 2.856954e-02        NA
# 11472 ENSG00000198682.12   PAPSS2 0.03138251 0.9540615 0.4061830 3.727603e-03        NA
# 11714 ENSG00000151892.14    GFRA1 0.04304279 0.9120729 0.8894487 4.078983e-02 0.8132686
# 12375  ENSG00000126500.3    FLRT1 0.03138251 0.8959630 0.4866056 6.467140e-06 0.5278196
# 16078 ENSG00000103653.16      CSK 0.04304279 0.1927468 0.3840742 4.854267e-03 0.8719310
# 17765 ENSG00000072310.16   SREBF1 0.03138251 0.6662341 0.9649995 1.223075e-02 0.5370652
# 18290 ENSG00000173868.11 PHOSPHO1 0.04304279 0.6822733 0.8546697 3.060767e-06        NA
# 22010 ENSG00000131831.17     RAI2 0.03138251 0.9752209 0.8936360 3.721027e-04 0.9014092


addmargins(table(all_FDR$Hb < 0.05, all_FDR$DG < 0.05, useNA = "ifany"))
#       FALSE  TRUE  <NA>   Sum
# FALSE 17867    10  4834 22711
# TRUE     25     0    20    45
# Sum   17892    10  4854 22756

pdf(file.path(dir_plots, "ggpairs_FDR.pdf"), height = 10, width = 10)
ggpairs(
    mutate(all_FDR,
        adj.P.Val = Hb,
        Hb = -log10(Hb),
        DLPFC = -log10(DLPFC),
        HIPPO = -log10(HIPPO),
        Caudate = -log10(Caudate),
        DG = -log10(DG)
    ),
    columns = c("Hb", "DLPFC", "HIPPO", "Caudate", "DG"),
    ggplot2::aes(
        colour = adj.P.Val < 0.05,
        alpha = ifelse(adj.P.Val < 0.05, 1, 1 / 3)
    )
) + theme_bw()
dev.off()


## Merge data: t-statistics
all_t <- mutate(habenula, Hb = t)[, c("gencodeID", "Symbol", "Hb", "adj.P.Val")]
all_t <- left_join(all_t, mutate(bsp2_dlpfc, DLPFC = t)[, c("gencodeID", "DLPFC")], by = "gencodeID")
all_t <- left_join(all_t, mutate(bsp2_hpc, HIPPO = t)[, c("gencodeID", "HIPPO")], by = "gencodeID")
all_t <- left_join(all_t, mutate(bsp3_caudate, Caudate = t)[, c("gencodeID", "Caudate")], by = "gencodeID")
all_t <- left_join(all_t, mutate(dg, DG = SZ_t)[, c("gencodeID", "DG")], by = "gencodeID")

p <- ggpairs(
    all_t,
    columns = c("Hb", "DLPFC", "HIPPO", "Caudate", "DG"),
    ggplot2::aes(
        colour = adj.P.Val < 0.05,
        alpha = ifelse(adj.P.Val < 0.05, 1, 1 / 3)
    )
) + theme_bw()
pdf(file.path(dir_plots, "ggpairs_t-stats.pdf"), height = 10, width = 10)
print(p)
dev.off()


set.seed(20231220)
p1 <- getPlot(p, 2, 1) +
    xlab("Hb SCZD vs control t-stat") +
    ylab("DLPFC SCZD vs control t-stat") +
    theme_bw(base_size = 20) +
    theme(legend.position = "none") +
    geom_text_repel(
        aes(label = Symbol),
        data = all_t[which(all_FDR$Hb < 0.05 & all_FDR$DLPFC < 0.05), ]
    )
p2 <- getPlot(p, 3, 1) +
    xlab("Hb SCZD vs control t-stat") +
    ylab("HIPPO SCZD vs control t-stat") +
    theme_bw(base_size = 20) +
    theme(legend.position = "none") +
    geom_text_repel(
        aes(label = Symbol),
        data = all_t[which(all_FDR$Hb < 0.05 & all_FDR$HIPPO < 0.05), ]
    )
p3 <- getPlot(p, 4, 1) +
    xlab("Hb SCZD vs control t-stat") +
    ylab("Caudate SCZD vs control t-stat") +
    theme_bw(base_size = 20) +
    theme(legend.position = "none") +
    geom_text_repel(
        aes(label = Symbol),
        data = all_t[which(all_FDR$Hb < 0.05 & all_FDR$Caudate < 0.05), ],
        force = 5
    )
p4 <- getPlot(p, 5, 1) +
    xlab("Hb SCZD vs control t-stat") +
    ylab("DG SCZD vs control t-stat") +
    theme_bw(base_size = 20) +
    scale_alpha(guide = 'none') +
    theme(
        legend.position = c(0.89, 0.11),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.background = element_rect(fill=alpha('grey', 0.2))
    ) +
    guides(colour = guide_legend(title =
            "Hb FDR < 0.05"))

pdf(
    file.path(dir_plots, "ggpairs_t-stats_Hb.pdf"),
    height = 10,
    width = 10
)
p1 + p2 + p3 + p4
dev.off()

pdf(
    file.path(dir_plots, "ggpairs_t-stats_Hb_facet.pdf"),
    height = 10,
    width = 10
)
p1 +
    facet_grid(reformulate('"Hb"', '"DLPFC"'), switch = "y") +
    xlab("") +
    ylab("") +
p2 +
    facet_grid(reformulate('"Hb"', '"HIPPO"')) +
    xlab("") +
    ylab("") +
p3 +
    facet_grid(reformulate('"Hb"', '"Caudate"'), switch = "both") +
    xlab("") +
    ylab("") +
p4 +
    facet_grid(reformulate('"Hb"', '"DG"'), switch = "x") +
    xlab("") +
    ylab("")
dev.off()

pdf(file.path(dir_plots, "cor_Hb_sig_vs_Hb_notSig_t-stats.pdf"))
plot_cor_type(all_t)
dev.off()
# 	One Sample t-test
#
# data:  df[, c("Sig", "NotSig")]
# t = 3.3082, df = 7, p-value = 0.01297
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#  0.07492158 0.45043804
# sample estimates:
# mean of x
# 0.2626798


## Merge data: log FC
all_logFC <- mutate(habenula, Hb = logFC)[, c("gencodeID", "Symbol", "Hb", "adj.P.Val")]
all_logFC <- left_join(all_logFC, mutate(bsp2_dlpfc, DLPFC = logFC)[, c("gencodeID", "DLPFC")], by = "gencodeID")
all_logFC <- left_join(all_logFC, mutate(bsp2_hpc, HIPPO = logFC)[, c("gencodeID", "HIPPO")], by = "gencodeID")
all_logFC <- left_join(all_logFC, mutate(bsp3_caudate, Caudate = logFC)[, c("gencodeID", "Caudate")], by = "gencodeID")
all_logFC <- left_join(all_logFC, mutate(dg, DG = SZ_logFC)[, c("gencodeID", "DG")], by = "gencodeID")

pdf(file.path(dir_plots, "ggpairs_logFC.pdf"), height = 10, width = 10)
ggpairs(
    all_logFC,
    columns = c("Hb", "DLPFC", "HIPPO", "Caudate", "DG"),
    ggplot2::aes(
        colour = adj.P.Val < 0.05,
        alpha = ifelse(adj.P.Val < 0.05, 1, 1 / 3)
    )
) + theme_bw()
dev.off()

pdf(file.path(dir_plots, "cor_Hb_sig_vs_Hb_notSig_logFC.pdf"))
plot_cor_type(all_logFC)
dev.off()
# 	One Sample t-test
#
# data:  df[, c("Sig", "NotSig")]
# t = 4.024, df = 7, p-value = 0.005034
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#  0.1086035 0.4181199
# sample estimates:
# mean of x
# 0.2633617

#### plot with overlap colors ####
library(tidyverse)

all_FDR_long <- all_FDR |>
    pivot_longer(!c(gencodeID:Hb), names_to = "region", values_to = "region_FDR") |>
    rename(Hb_FDR = Hb) |>
    mutate(signif = case_when(Hb_FDR < 0.05 & region_FDR < 0.05 ~ "Both",
                              region_FDR < 0.05 ~ "only_region",
                              Hb_FDR < 0.05 ~ "only_Hb",
                              TRUE ~ "None"))

all_FDR_long |> count(region, signif)

all_t_long <- all_t |>
    pivot_longer(!c(gencodeID:adj.P.Val), names_to = "region", values_to = "region_t") |>
    rename(Hb_t = Hb) |>
    left_join(all_FDR_long) |>
    mutate(text = signif == "Both" | Symbol %in% c("CCDC141", "QPRT", "HES5", "EHMT2"))

all_t_long |> filter(text) |> count(Symbol, signif)

tstat_scatter <-l

# ggsave(tstat_scatter, filename = file.path(dir_plots, "Hb_v_region_t-stats_scater_facet.pdf"), height = 5, width = 7)
ggsave(tstat_scatter, filename = file.path(dir_plots, "Hb_v_region_t-stats_scater_facet.png"), height = 5, width = 6)

## Save for later
save(
    all_FDR, all_t, all_logFC,
    file = file.path(dir_rdata, "SCZD_vs_control_Hb_and_other_regions.Rdata")
)

extract_info <- function(input, suffix) {
    stopifnot(identical(habenula$gencodeID, input$gencodeID))
    x <- input[, c("DLPFC", "HIPPO", "Caudate", "DG")]
    colnames(x) <- paste0(colnames(x), "_", suffix)
    return(x)
}
stopifnot(identical(habenula$Symbol, habenula$MGI_Symbol))
habenula_combined <- rename(
    cbind(
        habenula[, -which(
            colnames(habenula) %in% c(
                "r.name",
                "seqnames",
                "start",
                "end",
                "width",
                "strand",
                "NumTx",
                "gencodeTx",
                "Class",
                "MGI_Symbol",
                "EntrezID",
                "ensemblID"
            )
        )],
        extract_info(all_t, "t"),
        extract_info(all_logFC, "logFC"),
        extract_info(all_FDR, "adj.P.Val")
    ),
    Hb_meanExprs = meanExprs,
    Hb_logFC = logFC,
    Hb_AveExpr = AveExpr,
    Hb_t = t,
    Hb_P.Value = P.Value,
    Hb_adj.P.Val = adj.P.Val,
    Hb_B = B
)
write.csv(
    habenula_combined,
    row.names = FALSE,
    quote = FALSE,
    file = file.path(dir_rdata, "TableSxx_SZCD_vs_Control_gene_results.csv")
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ─────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.3.1 (2023-06-16)
#  os       macOS Sonoma 14.2
#  system   aarch64, darwin20
#  ui       RStudio
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       America/Mexico_City
#  date     2023-12-21
#  rstudio  2023.06.1+524 Mountain Hydrangea (desktop)
#  pandoc   3.1.5 @ /opt/homebrew/bin/pandoc
#
# ─ Packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package       * version date (UTC) lib source
#  BiocFileCache * 2.8.0   2023-04-25 [1] Bioconductor
#  bit             4.0.5   2022-11-15 [1] CRAN (R 4.3.0)
#  bit64           4.0.5   2020-08-30 [1] CRAN (R 4.3.0)
#  blob            1.2.4   2023-03-17 [1] CRAN (R 4.3.0)
#  brio            1.1.3   2021-11-30 [1] CRAN (R 4.3.0)
#  cachem          1.0.8   2023-05-01 [1] CRAN (R 4.3.0)
#  callr           3.7.3   2022-11-02 [1] CRAN (R 4.3.0)
#  cellranger      1.1.0   2016-07-27 [1] CRAN (R 4.3.0)
#  cli             3.6.1   2023-03-23 [1] CRAN (R 4.3.0)
#  colorout        1.3-0   2023-09-28 [1] Github (jalvesaq/colorout@8384882)
#  colorspace      2.1-0   2023-01-23 [1] CRAN (R 4.3.0)
#  crayon          1.5.2   2022-09-29 [1] CRAN (R 4.3.0)
#  curl            5.0.2   2023-08-14 [1] CRAN (R 4.3.1)
#  data.table      1.14.8  2023-02-17 [1] CRAN (R 4.3.0)
#  DBI             1.1.3   2022-06-18 [1] CRAN (R 4.3.0)
#  dbplyr        * 2.3.3   2023-07-07 [1] CRAN (R 4.3.0)
#  devtools      * 2.4.5   2022-10-11 [1] CRAN (R 4.3.0)
#  digest          0.6.33  2023-07-07 [1] CRAN (R 4.3.0)
#  dplyr         * 1.1.3   2023-09-03 [1] CRAN (R 4.3.0)
#  ellipsis        0.3.2   2021-04-29 [1] CRAN (R 4.3.0)
#  fansi           1.0.4   2023-01-22 [1] CRAN (R 4.3.0)
#  farver          2.1.1   2022-07-06 [1] CRAN (R 4.3.0)
#  fastmap         1.1.1   2023-02-24 [1] CRAN (R 4.3.0)
#  filelock        1.0.2   2018-10-05 [1] CRAN (R 4.3.0)
#  fs              1.6.3   2023-07-20 [1] CRAN (R 4.3.0)
#  generics        0.1.3   2022-07-05 [1] CRAN (R 4.3.0)
#  GGally        * 2.1.2   2021-06-21 [1] CRAN (R 4.3.0)
#  ggplot2       * 3.4.3   2023-08-14 [1] CRAN (R 4.3.0)
#  ggrepel       * 0.9.3   2023-02-03 [1] CRAN (R 4.3.0)
#  glue            1.6.2   2022-02-24 [1] CRAN (R 4.3.0)
#  gtable          0.3.4   2023-08-21 [1] CRAN (R 4.3.0)
#  here          * 1.0.1   2020-12-13 [1] CRAN (R 4.3.0)
#  hms             1.1.3   2023-03-21 [1] CRAN (R 4.3.0)
#  htmltools       0.5.6   2023-08-10 [1] CRAN (R 4.3.0)
#  htmlwidgets     1.6.2   2023-03-17 [1] CRAN (R 4.3.0)
#  httpuv          1.6.11  2023-05-11 [1] CRAN (R 4.3.0)
#  httr            1.4.7   2023-08-15 [1] CRAN (R 4.3.0)
#  labeling        0.4.3   2023-08-29 [1] CRAN (R 4.3.0)
#  later           1.3.1   2023-05-02 [1] CRAN (R 4.3.0)
#  lifecycle       1.0.3   2022-10-07 [1] CRAN (R 4.3.0)
#  lubridate       1.9.2   2023-02-10 [1] CRAN (R 4.3.0)
#  magrittr        2.0.3   2022-03-30 [1] CRAN (R 4.3.0)
#  memoise         2.0.1   2021-11-26 [1] CRAN (R 4.3.0)
#  mime            0.12    2021-09-28 [1] CRAN (R 4.3.0)
#  miniUI          0.1.1.1 2018-05-18 [1] CRAN (R 4.3.0)
#  munsell         0.5.0   2018-06-12 [1] CRAN (R 4.3.0)
#  patchwork     * 1.1.3   2023-08-14 [1] CRAN (R 4.3.1)
#  pillar          1.9.0   2023-03-22 [1] CRAN (R 4.3.0)
#  pkgbuild        1.4.2   2023-06-26 [1] CRAN (R 4.3.0)
#  pkgconfig       2.0.3   2019-09-22 [1] CRAN (R 4.3.0)
#  pkgload         1.3.2.1 2023-07-08 [1] CRAN (R 4.3.0)
#  plyr            1.8.8   2022-11-11 [1] CRAN (R 4.3.0)
#  prettyunits     1.1.1   2020-01-24 [1] CRAN (R 4.3.0)
#  processx        3.8.2   2023-06-30 [1] CRAN (R 4.3.0)
#  profvis         0.3.8   2023-05-02 [1] CRAN (R 4.3.0)
#  progress        1.2.2   2019-05-16 [1] CRAN (R 4.3.0)
#  promises        1.2.1   2023-08-10 [1] CRAN (R 4.3.0)
#  prompt          1.0.2   2023-08-31 [1] CRAN (R 4.3.0)
#  ps              1.7.5   2023-04-18 [1] CRAN (R 4.3.0)
#  purrr           1.0.2   2023-08-10 [1] CRAN (R 4.3.0)
#  R6              2.5.1   2021-08-19 [1] CRAN (R 4.3.0)
#  RColorBrewer    1.1-3   2022-04-03 [1] CRAN (R 4.3.0)
#  Rcpp            1.0.11  2023-07-06 [1] CRAN (R 4.3.0)
#  readxl        * 1.4.3   2023-07-06 [1] CRAN (R 4.3.0)
#  remotes         2.4.2.1 2023-07-18 [1] CRAN (R 4.3.0)
#  reshape         0.8.9   2022-04-12 [1] CRAN (R 4.3.0)
#  rlang           1.1.1   2023-04-28 [1] CRAN (R 4.3.0)
#  rprojroot       2.0.3   2022-04-02 [1] CRAN (R 4.3.0)
#  RSQLite         2.3.1   2023-04-03 [1] CRAN (R 4.3.0)
#  rsthemes        0.4.0   2023-05-06 [1] Github (gadenbuie/rsthemes@34a55a4)
#  rstudioapi      0.15.0  2023-07-07 [1] CRAN (R 4.3.0)
#  scales          1.2.1   2022-08-20 [1] CRAN (R 4.3.0)
#  sessioninfo   * 1.2.2   2021-12-06 [1] CRAN (R 4.3.0)
#  shiny           1.7.5   2023-08-12 [1] CRAN (R 4.3.0)
#  stringi         1.7.12  2023-01-11 [1] CRAN (R 4.3.0)
#  stringr         1.5.0   2022-12-02 [1] CRAN (R 4.3.0)
#  suncalc         0.5.1   2022-09-29 [1] CRAN (R 4.3.0)
#  testthat      * 3.1.10  2023-07-06 [1] CRAN (R 4.3.0)
#  tibble          3.2.1   2023-03-20 [1] CRAN (R 4.3.0)
#  tidyselect      1.2.0   2022-10-10 [1] CRAN (R 4.3.0)
#  timechange      0.2.0   2023-01-11 [1] CRAN (R 4.3.0)
#  urlchecker      1.0.1   2021-11-30 [1] CRAN (R 4.3.0)
#  usethis       * 2.2.2   2023-07-06 [1] CRAN (R 4.3.0)
#  utf8            1.2.3   2023-01-31 [1] CRAN (R 4.3.0)
#  vctrs           0.6.3   2023-06-14 [1] CRAN (R 4.3.0)
#  withr           2.5.0   2022-03-03 [1] CRAN (R 4.3.0)
#  xtable          1.8-4   2019-04-21 [1] CRAN (R 4.3.0)
#
#  [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
#
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
