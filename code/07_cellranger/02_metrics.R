metrics_files <-
    list.files(
        here::here("processed-data", "07_cellranger"),
        pattern = "metrics_summary.csv",
        full.names = TRUE,
        recursive = TRUE
    )
names(metrics_files) <-
    list.files(here::here("processed-data", "07_cellranger"), pattern = "^Br")

metrics <- do.call(rbind, lapply(metrics_files, read.csv))

sum(as.numeric(gsub(",", "", metrics$Estimated.Number.of.Cells)))
# 20327

summary(as.numeric(gsub(",", "", metrics$Estimated.Number.of.Cells)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  923    2024    3156    2904    3906    4389

summary(as.numeric(gsub(",", "", metrics$Number.of.Reads)) / 1e6)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 162.8   172.9   193.0   188.3   201.8   213.0

summary(as.numeric(gsub(",", "", metrics$Mean.Reads.per.Cell)))
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 39476   45998   54641   88719  102062  230795

summary(as.numeric(gsub(
    ",", "", metrics$Median.UMI.Counts.per.Cell
)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 7752   10856   14135   17674   22824   34469

summary(as.numeric(gsub(",", "", metrics$Median.Genes.per.Cell)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 2836    3838    4592    5083    6316    7843

write.csv(
    cbind(data.frame(BrNum = rownames(metrics)), metrics),
    file = here::here(
        "processed-data",
        "07_cellranger",
        "STable_CellRanger_metrics.csv"
    ),
    row.names = FALSE
)
