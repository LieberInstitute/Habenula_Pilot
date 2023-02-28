## 2/16/23 - Bukola Ajanaku
# Labeling clusters based on gene expression plots. This script loads the excel sheet where this 
# info is saved and then adds it to the official sce object so far (sce_post_clustering_with_logcounts),
# not to be confused with the pseudobulked sce objects as that is just for heatmap purposes :)
# qrsh -l mem_free=20G,h_vmem=20G
# Terminal 1

library("SingleCellExperiment")
library("here")
library("xlsx")
library("tidyverse")

# loading sce object with clustered harmonized data
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
          "sce_post_clustering_with_logcounts.Rdata"))

# pulling excel sheet with annotations
annoWT10 <- read.xlsx(file = here("processed-data", "04_snRNA-seq", "07b_Marking_Clusters",
                      "GeneAnnotations_Bukola.xlsx"), sheetName = c("WalkTrap10"))
annoWT20 <- read.xlsx(file = here("processed-data", "04_snRNA-seq", "07b_Marking_Clusters",
                                  "GeneAnnotations_Bukola.xlsx"), sheetName = c("WalkTrap20"))
annoWT50 <- read.xlsx(file = here("processed-data", "04_snRNA-seq", "07b_Marking_Clusters",
                                  "GeneAnnotations_Bukola.xlsx"), sheetName = c("WalkTrap 50"))

# cleaning up for standardization
# 10
annoWT10_clean <- annoWT10 |> 
  filter(!is.na(Cluster)) |>
  mutate(Type_clean = gsub("[^a-zA-Z0-9]", "", Type), 
         Cluster = paste0("10wTrap_", Cluster))

#20
annoWT20_clean <- annoWT20 |> 
  filter(!is.na(Cluster)) |>
  mutate(Type_clean = gsub("[^a-zA-Z0-9]", "", Type), 
         Cluster = paste0("20wTrap_", Cluster))

#30
annoWT50_clean <- annoWT50 |> 
  filter(!is.na(Cluster)) |>
  mutate(Type_clean = gsub("[^a-zA-Z0-9]", "", Type), 
         Cluster = paste0("50wTrap_", Cluster))

# sanity check (looking for no repeat names in summarized list)
annoWT10_clean |> count(Type_clean)
annoWT20_clean |> count(Type_clean)
annoWT50_clean |> count(Type_clean)

# using match to add annotated name columns for wt50
sce$cellType_wT10 <- annoWT10_clean$Type_clean[match(sce$wT_10_Erik, annoWT10_clean$Cluster)]
sce$cellType_wT20 <- annoWT20_clean$Type_clean[match(sce$wT_20_Erik, annoWT20_clean$Cluster)]
sce$cellType_wT50 <- annoWT50_clean$Type_clean[match(sce$wT_50_Erik, annoWT50_clean$Cluster)]

# saving

save(sce, file = here("processed-data", "04_snRNA-seq", "sce_objects", 
                      "sce_post_clustering_with_logcounts.Rdata"))
