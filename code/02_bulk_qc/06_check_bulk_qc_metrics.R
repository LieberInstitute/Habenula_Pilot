## Louise Huuki-Myers, Feb 2025
## Check quality metrics vs. other LIBD datasets

library("SummarizedExperiment")
library("tidyverse")
library("sessioninfo")
library("here")
library("readxl")
library("ggrepel")
library("jaffelab")
library("scuttle")
library("GGally")
library("patchwork")
library("HDF5Array")

## prep dirs ##
plot_dir <- here("plots", "02_bulk_qc", "06_check_bulk_qc_metrics")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load Big RSE Data ####
rse_gene <- loadHDF5SummarizedExperiment("/dcs04/lieber/lcolladotor/dbDev_LIBD001/RNAseq/datasets/brain_swap/dataset_RSEs/h5_rse_gene")
dim(rse_gene)
# [1] 58037  6243

## missing some details about library type
pd_big <- as.data.frame(colData(rse_gene)) |>
    mutate(library_type = ifelse(Dataset == "BrainSeq_Phase1", "polyA", "RiboZero"))

table(pd_big$library_type)
# polyA RiboZeroGold
# 746         5497

table(pd_big$Dataset)
# Astellas_DG   BrainSeq_Phase2_DLPFC   BrainSeq_Phase2_HIPPO BrainSeq_Phase3_Caudate     BrainSeq_Phase4and5
# 263                     453                     447                     464                     490
# Habenula            Nicotine_NAc        psychENCODE_Mood                 VA_PTSD         BrainSeq_Phase1
# 69                     235                    1108                    1280                     746
# PTSD_BrainOmics
# 688

colnames(pd_big)
# [1] "SAMPLE_ID"         "RNum"              "RIN"               "Region"            "Dataset"
# [6] "BrNum"             "Dx"                "Age"               "Sex"               "Race"
# [11] "numReads"          "numMapped"         "numUnmapped"       "mitoMapped"        "totalMapped"
# [16] "overallMapRate"    "concordMapRate"    "mitoRate"          "rRNA_rate"         "totalAssignedGene"
# [21] "bamFile"           "rna_preSwap_BrNum" "library_type"

## historic cutoffs

historic_cutoffs <- pd_big |>
    group_by(library_type, Dataset) |>
    summarise(
        concordMapRate_min = min(concordMapRate),
        mitoMapped_max = max(mitoMapped),
        mitoRate_max = max(mitoRate),
        numMapped_min = min(numMapped),
        numReads_min = min(numReads),
        overallMapRate_min = min(overallMapRate),
        rRNArate_max = max(rRNA_rate),
        totalAssignedGene_min = min(totalAssignedGene),
        totalMapped_min = min(totalMapped)
    )

historic_cutoff_long <- historic_cutoffs |>
    pivot_longer(!c(library_type, Dataset), names_to = "cutoff")

historic_cutoff_box <- historic_cutoff_long |>
    ggplot(aes(x = library_type, y = value)) +
    geom_boxplot(outlier.shape = NA) +
    geom_text_repel(aes(label = Dataset), color = "blue", size = 1) +
    facet_wrap(~cutoff, scales = "free_y") +
    theme_bw()

ggsave(historic_cutoff_box, filename = here(plot_dir, "historic_cutoff_boxplot.png"))

## sanity check numbers
historic_cutoff_median <- historic_cutoff_long |>
    filter(library_type == "RiboZero") |>
    group_by(library_type, cutoff) |>
    summarise(Median = median(value)) |>
    separate(cutoff, into = c("qc_var", "type"), sep = "_") |>
    mutate(qc_var = gsub("rRNArate", "rRNA_rate", qc_var))

# historic_cutoff_median <- read.csv(here("processed-data", "02_quality_control","QC_cutoffs_historic.csv"))

# qc_var            type   Median
# <chr>             <chr>   <dbl>
# 1 concordMapRate    min   5.22e-1
# 2 mitoMapped        max   1.19e+7
# 3 mitoRate          max   9.96e-2
# 4 numMapped         min   4.93e+7
# 5 numReads          min   6.75e+7
# 6 overallMapRate    min   5.47e-1
# 7 rRNA_rate         max   1.16e-5
# 8 totalAssignedGene min   2.94e-1
# 9 totalMapped       min   4.66e+7

#### Eval Hd data ####
pd_Hb <- pd_big |> filter(Dataset == "Habenula")
colnames(pd_Hb)

## Subset and QC data
qc_variables <- c("numReads", "numMapped", "numUnmapped", "overallMapRate", "concordMapRate", "totalMapped", "mitoMapped", "mitoRate", "rRNA_rate", "totalAssignedGene")

pd_qc_long <- pd_big |>
    select(SAMPLE_ID, Dataset,BrNum, library_type, all_of(qc_variables)) |>
    pivot_longer(!c(SAMPLE_ID, Dataset,BrNum, library_type), names_to = "qc_var") |>
    mutate(Habenula = Dataset == "Habenula") |>
    filter(library_type == "RiboZero")

pd_qc_long |> count(Habenula)

qc_boxplots_LIBD_exp <- ggplot(data = pd_qc_long, aes(x = Dataset, y = value)) +
    geom_boxplot(aes(fill = Habenula)) +
    facet_wrap(~qc_var, scales = "free_y", ncol = 5)  +
    geom_hline(data = historic_cutoff_median, aes(yintercept = Median), color = "darkgreen") +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 17)
    )

ggsave(qc_boxplots_LIBD_exp, filename = here(plot_dir, "qc_boxplots_LIBD_experiment_Hb.png"), width = 24, height = 12)

## overallMapRate ONLY
mapRate_boxplots_LIBD_exp <- ggplot(data = pd_qc_long |> filter(qc_var == "overallMapRate"),
                               aes(x = Dataset, y = value)) +
    geom_boxplot(aes(fill = Habenula)) +
    geom_text_repel(aes(label = ifelse(Habenula, BrNum, "")),
                     size = 1.7) +
    # geom_hline(data = historic_cutoff_median |> filter(qc_var == "overallMapRate"),
    #            aes(yintercept = Median), color = "darkgreen") +
    labs(y = "overallMapRate") +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 17)
    )

ggsave(mapRate_boxplots_LIBD_exp, filename = here(plot_dir, "mapRate_boxplots_LIBD_experiment_Hb.png"))


#### gg pairs ####
## how do the variables relate
focused_qc_metrics <- c(
    "concordMapRate",
    "mitoRate",
    "numMapped",
    "numReads",
    "overallMapRate",
    "totalAssignedGene",
    "totalMapped"
)

pd_Hb |> select(all_of(focused_qc_metrics)) |> summary()
# concordMapRate      mitoRate         numMapped            numReads         overallMapRate   totalAssignedGene
# Min.   :0.5216   Min.   :0.01315   Min.   : 23775158   Min.   : 30136424   Min.   :0.5471   Min.   :0.3226
# 1st Qu.:0.7691   1st Qu.:0.02879   1st Qu.: 50344020   1st Qu.: 60119580   1st Qu.:0.8177   1st Qu.:0.4403
# Median :0.8043   Median :0.03319   Median : 77489327   Median : 90648196   Median :0.8506   Median :0.4626
# Mean   :0.7918   Mean   :0.03704   Mean   :128016012   Mean   :150071886   Mean   :0.8368   Mean   :0.4757
# 3rd Qu.:0.8454   3rd Qu.:0.04189   3rd Qu.:197911435   3rd Qu.:235646592   3rd Qu.:0.8938   3rd Qu.:0.5215
# Max.   :0.8716   Max.   :0.07412   Max.   :589967944   Max.   :645798490   Max.   :0.9169   Max.   :0.6105
# totalMapped
# Min.   : 24300119
# 1st Qu.: 52818444
# Median : 80534264
# Mean   :132484278
# 3rd Qu.:205381009
# Max.   :612178205




qc_ggpair <- ggpairs(pd_Hb,
                     columns = focused_qc_metrics) +
    theme_bw() +
    labs(title = "Habenula Bulk")

ggsave(qc_ggpair, filename = here(plot_dir, paste0("Hb_qc_ggpairs_focus.png")), width = 10, height = 10)


# #### Combined Data automatic outliers ####
#
# names(focused_qc_metrics) <- focused_qc_metrics
# tail <- c(
#     concordMapRate = "lower",
#     mitoRate = "higher",
#     numMapped = "lower",
#     numReads = "lower",
#     overallMapRate = "lower",
#     totalAssignedGene = "lower",
#     totalMapped = "lower",
#     mitoRate = "higher"
# )
#
# tail <- tail[names(focused_qc_metrics)]
#
# nmads <- 3
# qc_guide <- tibble(
#     qc_var = focused_qc_metrics,
#     tail = tail,
#     nmads = nmads
# )
#
# auto_outliers <- pmap(qc_guide, function(...) {
#     r <- tibble(...)
#     out <- isOutlier(rse_gene[[r$qc_var]], batch = rse_gene$library_type, type = r$tail, nmads = r$nmads)
#     return(out)
# })
#
#
# ## find MAD cutoffs
# auto_cutoff <- do.call(
#     "rbind",
#     map2(auto_outliers, names(auto_outliers), ~ attr(.x, "thresholds") |>
#              as.data.frame() |>
#              rownames_to_column("cutoff_type") |>
#              add_column(qc_var = .y, .before = 1))
# ) |>
#     rename(cutoff = 3)
#
# auto_cutoff |>
#     filter(abs(cutoff) != Inf)
#
# #                                qc_var cutoff_type       cutoff
# # concordMapRate.1       concordMapRate       lower 5.737424e-01
# # mitoRate.2                   mitoRate      higher 6.776782e-02
# # numMapped.1                 numMapped       lower 3.348988e+07
# # numReads.1                   numReads       lower 5.375371e+07
# # overallMapRate.1       overallMapRate       lower 6.531191e-01
# # totalAssignedGene.1 totalAssignedGene       lower 1.811326e-01
# # totalMapped.1             totalMapped       lower 2.992672e+07
#
#
# map(auto_outliers, table)
#
# drop_metrics <- c("numMapped", "numReads", "totalMapped")
#
# auto_outliers_drop <- Reduce("|", auto_outliers[drop_metrics]) ## clearly poor metrics, also catches extremes for overallMapRate & concordMapRate
# table(auto_outliers_drop)
# # FALSE  TRUE
# # 111     2
#
# auto_outliers_warn <- auto_outliers$totalAssignedGene | (auto_outliers$mitoRate & rse_gene$library_type == "polyA") ## mitoRate for polyA, totalAssignedGene for RiboZero?
#
# auto_outlier_tb <- tibble(
#     SAMPLE_ID = pd_Hb$SAMPLE_ID,
#     auto_drop = auto_outliers_drop,
#     auto_warn = auto_outliers_warn
# ) |>
#     mutate(qc_class = ifelse(auto_drop, "drop", ifelse(auto_warn, "warn", "pass")))
#
# ## save in csv file
# write.csv(auto_outlier_tb,
#           file = here("processed-data", "02_quality_control", "QC_record_DLPFC_bulk.csv"),
#           row.names = FALSE
# )
#
# auto_outlier_tb |>
#     filter(auto_drop | auto_warn)
#
# ## add to new qc_long
# pd_Hb_qc_long <- pd_Hb_qc_long |>
#     left_join(auto_outlier_tb)
#
# ## selected cutoffs
# auto_cutoff_select <- auto_cutoff |>
#     mutate(qc_class = ifelse(qc_var %in% drop_metrics, "drop",
#                              ifelse(qc_var %in% c("mitoRate", "totalAssignedGene"), "warn", NA)
#     )) |>
#     filter(
#         abs(cutoff) != Inf,
#         !is.na(qc_class),
#         (qc_var != "mitoRate" | library_type != "RiboZeroGold")
#     ) |>
#     arrange(qc_class)
#
# # qc_var            cutoff_type library_type       cutoff qc_class
# # <chr>             <chr>       <chr>               <dbl> <chr>
# #   1 numMapped         lower       polyA        61280789.    drop
# # 2 numMapped         lower       RiboZeroGold 48838795.    drop
# # 3 numReads          lower       polyA        65036118.    drop
# # 4 numReads          lower       RiboZeroGold 51168785.    drop
# # 5 totalMapped       lower       polyA        45653338.    drop
# # 6 totalMapped       lower       RiboZeroGold 50556432.    drop
# # 7 mitoRate          higher      polyA               0.340 warn
# # 8 totalAssignedGene lower       polyA               0.548 warn
# # 9 totalAssignedGene lower       RiboZeroGold        0.166 warn
#
# ## record cutoffs
# write.csv(auto_cutoff_select,
#           file = here("processed-data", "02_quality_control", "QC_cutoffs_DLPFC_bulk.csv"),
#           row.names = FALSE
# )
#
# ### boxplot for global isOutlier cutoffs
# qc_boxplot_isOutlier <- pd_Hb_qc_long |>
#     filter(qc_var %in% focused_qc_metrics) |>
#     ggplot(aes(x = library_type, y = value)) +
#     geom_boxplot(aes(fill = library_type), outlier.shape = NA, alpha = 0.5) +
#     geom_jitter(aes(shape = qc_class)) +
#     facet_wrap(~qc_var, scales = "free_y", nrow = 2) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     labs(title = paste("isOutlier cutoff -", nmads, "MADs")) +
#     # scale_color_manual(values = c(`FALSE` = "black", `TRUE` = "grey50")) +
#     # geom_hline(data = auto_cutoff |> filter(abs(cutoff) != Inf),
#     #            aes(yintercept = cutoff, color = library_type), linetype = "dashed") +
#     geom_hline(
#         data = auto_cutoff_select,
#         aes(yintercept = cutoff, color = library_type, linetype = qc_class)
#     ) +
#     geom_hline(
#         data = historic_cutoff_median |> filter(qc_var %in% focused_qc_metrics),
#         aes(yintercept = Median), color = "darkgreen"
#     )
# #
# ggsave(qc_boxplot_isOutlier, filename = here(plot_dir, paste0("qc_boxplots_isOutlier.png")), width = 12)
#
#
# #### New Data Only Boxplots ####
#
# qc_boxplots_lt <- pd_Hb_qc_long |>
#     group_by(library_type) |>
#     group_map(
#         ~ {
#             qc_boxplot <- ggplot(.x, aes(x = Dataset, y = value)) +
#                 geom_boxplot(outlier.shape = NULL) +
#                 geom_jitter(aes(shapre = qc_class), width = .2) +
#                 facet_wrap(~qc_var, scales = "free_y", nrow = 2) +
#                 scale_color_manual(values = c(new = "black", old = "grey50")) +
#                 theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#                 labs(title = .y) +
#                 geom_hline(
#                     data = historic_cutoff_median |> filter(qc_var %in% focused_qc_metrics),
#                     aes(yintercept = Median), color = "darkgreen"
#                 )
#
#             ggsave(qc_boxplot, filename = here(plot_dir, paste0("qc_boxplots_", .y, ".png")), width = 12)
#             return(qc_boxplot)
#         }
#     )
#
# ggsave(qc_boxplots_lt[[1]] / qc_boxplots_lt[[2]], filename = here(plot_dir, "qc_boxplots.png"), width = 24, height = 12)
#
# ## color by library_prep
# pd_Hb_qc_long |>
#     group_by(library_type) |>
#     group_map(
#         ~ {
#             qc_boxplot <- ggplot(.x, aes(x = Dataset, y = value)) +
#                 geom_boxplot(aes(fill = library_prep), outlier.shape = NA, alpha = 0.5) +
#                 geom_point(aes(color = library_prep, shape = qc_class), position = position_jitterdodge()) +
#                 geom_text(aes(label = ifelse(qc_class == "drop", Sample, "")), size = 2) +
#                 facet_wrap(~qc_var, scales = "free_y", nrow = 2) +
#                 theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#                 labs(title = .y) +
#                 geom_hline(
#                     data = auto_cutoff_select |>
#                         filter(library_type == as.character(.y)),
#                     aes(yintercept = cutoff, linetype = qc_class)
#                 )
#
#             ggsave(qc_boxplot, filename = here(plot_dir, paste0("qc_boxplots_", .y, ".png")), width = 12)
#             # return(qc_boxplot)
#         }
#     )
#
# #### ERCC data ####
# pd_Hb <- pd_Hb |>
#     left_join(auto_outlier_tb)
#
# ercc_boxplot <- pd_Hb |>
#     ggplot(aes(x = Dataset, y = ERCCsumLogErr, fill = library_type)) +
#     geom_boxplot(outlier.shape = NA, alpha = 0.5) +
#     geom_jitter(aes(color = qc_class), width = .2) +
#     scale_color_manual(values = c(`FALSE` = "black", `TRUE` = "red")) +
#     facet_wrap(~round, scales = "free") +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     labs(subtitle = "No ERCC spike in for Round 1")
#
# ggsave(ercc_boxplot, filename = here(plot_dir, "ERCCsumLogErr_boxplot.png"))
#
# ## there was no ERCC spike in for round1 so those values are NA
#
# ## trend over library_prep
# ercc_boxplot_prep <- pd_Hb |>
#     filter(round == 2) |>
#     ggplot(aes(x = library_prep, y = ERCCsumLogErr, fill = library_prep)) +
#     geom_boxplot(outlier.shape = NA, alpha = 0.5) +
#     geom_jitter(aes(color = qc_class), width = .2) +
#     scale_color_manual(values = c(`FALSE` = "black", `TRUE` = "red")) +
#     facet_wrap(~library_type)
#
# ggsave(ercc_boxplot_prep, filename = here(plot_dir, "ERCCsumLogErr_boxplot_prep.png"))
#
#
# #### Metrics vs. ERCC ###
# ercc_check <- pd_Hb |>
#     select(SAMPLE_ID, ERCCsumLogErr, library_prep, round, qc_class) |>
#     # filter(!is.na(ERCCsumLogErr)) |>
#     replace_na(list(ERCCsumLogErr = Inf)) |>
#     left_join(
#         pd_Hb_qc_long |>
#             filter(qc_var %in% focused_qc_metrics),
#         multiple = "all"
#     )
#
#
# ercc_check |>
#     group_by(library_type) |>
#     group_map(~ {
#         ercc_scatter <- ggplot(.x, aes(x = ERCCsumLogErr, y = value, color = library_prep)) +
#             geom_point(aes(shape = qc_class)) +
#             # scale_color_manual(values = c(`FALSE` = "black", `TRUE` = "red"))+
#             facet_wrap(~qc_var, scales = "free_y", nrow = 1) +
#             labs(title = .y) +
#             geom_hline(
#                 data = auto_cutoff_select |>
#                     filter(library_type == as.character(.y)),
#                 aes(yintercept = cutoff, linetype = qc_class)
#             )
#
#         ggsave(ercc_scatter, filename = here(plot_dir, paste0("ERCC_scatter_", .y, ".png")), width = 14, height = 5)
#     })


# sgejobs::job_single('06_check_bulk_qc_metrics', create_shell = TRUE, queue= 'bluejay', memory = '5G', command = "Rscript 06_check_bulk_qc_metrics.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
