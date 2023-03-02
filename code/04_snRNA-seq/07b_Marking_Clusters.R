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
library("fasthplus")


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

##### Runing fast h plus on cluster methods ####
## gotten from DLPFC code by Louise Huuki-Myers:
# (https://app.slack.com/client/T029RHR07/D044U4QSPAM/thread/G019V8X9CVC-1677597764.265859#:~:text=https%3A//github.com/LieberInstitute/DLPFC_snRNAseq/blob/cbd94e541678f715cbeda9f045df5adbd2c0ef96/code/03_build_sce/cluster_mb_kmeans.R%23L33%2DL59)
find_t <- function(L, proportion = 0.05) {
  initial_t <- floor(length(L) * proportion)
  # message("init. t = ", initial_t)
  smallest_cluster_size <- min(table(L))
  n_labels <- length(unique(L))
  message("smallest cluster: ", smallest_cluster_size, ", n lables: ", n_labels)
  ifelse(smallest_cluster_size > (initial_t / n_labels), initial_t, smallest_cluster_size * n_labels)
}

# Running per cluster  
message("Find fasthplus for clusters - ", Sys.time())

fasthplus <- lapply(c(sce$wT_10_Erik, sce$wT_20_Erik, sce$wT_50_Erik), function(li) {
  message(Sys.time())
  initial_t <- find_t(L = li, proportion = 0.01)
  h <- hpb(
    D = reducedDims(sce)$HARMONY,
    L = li,
    t = initial_t,
    r = 30
  )
})

fasthplus <- unlist(fasthplus)

# # saving
# save(sce, file = here("processed-data", "04_snRNA-seq", "sce_objects", 
#                       "sce_post_clustering_with_logcounts.Rdata"))
