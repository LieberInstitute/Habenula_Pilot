## 2/16/23 - Bukola Ajanaku
# Labeling clusters based on gene expression plots.
# qrsh -l mem_free=20G,h_vmem=20G

library("SingleCellExperiment")
library("here")
library("xlsx")
library("tidyverse")

# loading sce object with clustered harmonized data
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
          "sce_post_clustering_with_logcounts.Rdata"))

# pulling excel sheet with annotations
annoWT10 <- read.xlsx(file = here("processed-data", "04_snRNA-seq", "07b_Marking_Clusters",
                      "GeneAnnotation_Bukola.xlsx"), sheetName = c("Walktrap10"))
annoWT20 <- read.xlsx(file = here("processed-data", "04_snRNA-seq", "07b_Marking_Clusters",
                                  "GeneAnnotation_Bukola.xlsx"), sheetName = c("WalkTrap20"))
annoWT50 <- read.xlsx(file = here("processed-data", "04_snRNA-seq", "07b_Marking_Clusters",
                                  "GeneAnnotation_Bukola.xlsx"), sheetName = c("WalkTrap 50"))

# cleaning up annoWT50 for standardization
annoWT50_clean <- annoWT50 |> 
  filter(!is.na(Cluster)) |>
  mutate(Type_clean = gsub("[^a-zA-Z0-9]", "", Type), 
         Cluster = paste0("50wTrap_", Cluster))

# sanity check 
annoWT50_clean |> count(Type_clean)

# using match to add annotated name columns for wt50
sce$cellType_wT50 <- annoWT50_clean$Type_clean[match(sce$wT_50_Erik, annoWT50_clean$Cluster)]



# saving

save(sce, file = here("processed-data", "04_snRNA-seq", "sce_objects", 
                      "sce_post_clustering_with_logcounts.Rdata"))



# 