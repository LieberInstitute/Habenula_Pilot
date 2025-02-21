## Louise Huuki-Myers, Feb, 2025
## Check cellRanger metrics vs. other LIBD studies, part of review

library(tidyverse)
library(here)
library(ggrepel)

## prep dirs ##
plot_dir <- here("plots", "04_snRNA-seq", "10_libd_cellRanger_metrics")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### other snRNA-seq stats ####

## Tran 10x & spatialDLPFC (we already comapred metrics there)
tran_metrics <- read_csv("/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/02_cellranger_metrics/tran_metrics.csv")|>
    mutate(sample_id  = as.character(`...1`),
           dataset = "Tran10x") |>
    select(-`...1`)

spatialDLPFC_metrics <- read_csv("/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/02_cellranger_metrics/cellranger_metrics.csv") |>
    mutate(sample_id  = `...1`,
           dataset = "SpatialDLPFC")|>
    select(-`...1`)

setdiff(colnames(tran_metrics), colnames(spatialDLPFC_metrics))
setdiff(colnames(spatialDLPFC_metrics), colnames(tran_metrics))

## Hb
Hb_metrics <- read_csv(here("processed-data", "07_cellranger", "STable_CellRanger_metrics.csv")) |>
    mutate(sample_id  = BrNum,
           dataset = "Hb_pilot",
           across(Valid.Barcodes:Fraction.Reads.in.Cells,~as.numeric(gsub("%", "", .x)))
    )

colnames(Hb_metrics)

Hb_metrics |> select(BrNum, Valid.Barcodes, Sequencing.Saturation, Q30.Bases.in.Barcode, Q30.Bases.in.RNA.Read)

Hb_cols <- colnames(Hb_metrics)

bind_rows()

## All metrics

all_metrics <- Hb_metrics |>
    bind_rows(tran_metrics)|>
    bind_rows(spatialDLPFC_metrics)

all_metrics |> group_by(dataset) |> summarize(Mean.Reads.per.Cell = max(Mean.Reads.per.Cell))

colnames(all_metrics)

## Concern metrics and Samples
# Br5558 demonstrates low cell recovery (10%), as well as concerningly high Mean.Reads.per.Cell, high Median.Genes.per.Cell, and high Median.UMI.Counts.per.Cell
# Br1469 is also concerning for high Sequencing.Saturation and low Median.UMI.Counts.per.Cell

concern_metrics <- c("Mean.Reads.per.Cell", "Median.Genes.per.Cell", "Median.UMI.Counts.per.Cell","Sequencing.Saturation")
concern_samples <- c("Br5558", "Br1469")

qc_long <- all_metrics |> select(all_of(c("sample_id","dataset", concern_metrics))) |>
    pivot_longer(!c("sample_id","dataset"), names_to = "qc_metric", values_to = "qc_value") |>
    mutate(concern_sample = sample_id %in% concern_samples)

quality_boxplot <- qc_long |>
    ggplot(aes(x = dataset, y = qc_value, color = dataset)) +
    geom_boxplot(outlier.shape = NULL) +
    geom_jitter() +
    geom_text(aes(label = ifelse(concern_sample, sample_id, "")),
                    color = "black", size = 2) +
    theme_bw() +
    facet_wrap(~qc_metric, scales = "free_y") +
    theme(legend.position = "None") +
    labs( y = "Cell Ranger Metric Value")

ggsave(quality_boxplot, filename = here(plot_dir, "quality_boxplot.png"))


