library("here")
library("BiocFileCache")
library("readxl")
library("dplyr")
library("GGally")
library("sessioninfo")

## Output directories
dir_rdata <- here("processed-data", "14_compare_vs_BrainSeq")
dir_plots <- here("plots", "14_compare_vs_BrainSeq")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

## Inputs
dir_input <- here("processed-data", "10_DEA", "04_DEA")
# dir_rawdata <- "/dcs04/lieber/lcolladotor/qSVA_LIBD3080/qsva_brain/brainseq_phase2_qsv/rdas"
dir_rawdata <- here("raw-data", "14_compare_vs_BrainSeq")

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
    print(t.test(df[, c("Sig", "NotSig")]))
    ggplot(df, aes(x = NotSig, y = Sig, colour = Col)) +
        geom_point(size = 5) +
        xlab("cor with Hb FDR >= 0.05") +
        ylab("cor with Hb FDR < 0.05") +
        scale_color_manual(values = region_colors) +
        guides(colour = guide_legend(title="Region")) +
        geom_abline(intercept = 0, slope = 1, col = "red") +
        xlim(range(c(df$Sig, df$NotSig))) +
        ylim(range(c(df$Sig, df$NotSig))) +
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
)


## Merge data: t-statistics
all_t <- mutate(habenula, Hb = t)[, c("gencodeID", "Symbol", "Hb", "adj.P.Val")]
all_t <- left_join(all_t, mutate(bsp2_dlpfc, DLPFC = t)[, c("gencodeID", "DLPFC")], by = "gencodeID")
all_t <- left_join(all_t, mutate(bsp2_hpc, HIPPO = t)[, c("gencodeID", "HIPPO")], by = "gencodeID")
all_t <- left_join(all_t, mutate(bsp3_caudate, Caudate = t)[, c("gencodeID", "Caudate")], by = "gencodeID")
all_t <- left_join(all_t, mutate(dg, DG = SZ_t)[, c("gencodeID", "DG")], by = "gencodeID")

ggpairs(
    all_t,
    columns = c("Hb", "DLPFC", "HIPPO", "Caudate", "DG"),
    ggplot2::aes(
        colour = adj.P.Val < 0.05,
        alpha = ifelse(adj.P.Val < 0.05, 1, 1 / 3)
    )
)

plot_cor_type(all_t)
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

ggpairs(
    all_logFC,
    columns = c("Hb", "DLPFC", "HIPPO", "Caudate", "DG"),
    ggplot2::aes(
        colour = adj.P.Val < 0.05,
        alpha = ifelse(adj.P.Val < 0.05, 1, 1 / 3)
    )
)

plot_cor_type(all_logFC)
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


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
