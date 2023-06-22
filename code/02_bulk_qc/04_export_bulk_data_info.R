
library("SummarizedExperiment")
library("dplyr")
library("here")

load(here("processed-data","rse_objects","rse_gene_filt_DEA_n69.rda"), verbose = TRUE)
dim(rse_gene)

## drop bad deconvolution sample
"Br5572" %in% pd$BrNum
rse_gene <- rse_gene[,rse_gene$BrNum != "Br5572"]

pd <- as.data.frame(colData(rse_gene)) |>
  select(BrNum, RNum, PrimaryDx, AgeDeath, Race, Sex, Brain.Region)

write.csv(pd, here("processed-data", "02_bulk_qc","Habenula_pd_export_n68.csv"), row.names = FALSE)

pd |> count(PrimaryDx)

pd |>
  filter(PrimaryDx == "Control") |> select(BrNum) |>
  write.csv(here("processed-data", "02_bulk_qc", "Habenula_BrNum_control.csv"), row.names = FALSE)





