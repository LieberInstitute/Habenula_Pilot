## March 2, 2023 - Bukola Ajanaku
# # We've decided (based on heatmaps against Erik's clusters and against our gene
# # marker interest list) to proceed with Walktrap 10 which has 37  clusters.
# # This script will updated the annotations for each cluster based on further investigation
# # and also plot sn-RNA-seq qc against for each cluster to make sure no cluster is
# # driven by poor QC.
# qrsh -l mem_free=20G,h_vmem=20G

library("SingleCellExperiment")
library("here")
library("sessioninfo")
library("xlsx")
library("tidyverse")


# most updaated sce object with first round of annotations (where we selected to move forward with wT10)
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
          "sce_post_clustering_with_logcounts.Rdata"))

# pulling excel sheet with updated annotations (looking for snType and bulkType)
snAnno <- read.xlsx(file = here("processed-data", "04_snRNA-seq", "09_Clustered_QC",
                                  "GeneAnnotations_Bukola.xlsx"), sheetName = c("WalkTrap10"))

# cleaning up for standardization
snAnno_clean <- snAnno |> 
  filter(!is.na(Cluster)) |>
  mutate(Type_clean = gsub("[^a-zA-Z0-9]", "", Type), 
         Cluster = paste0("10wTrap_", Cluster)) 

# dropping three extra columns 
snAnno_clean <- select(snAnno_clean, -c(starts_with("Na"), "MoreInfo"))

# sanity check (checking for repeats in annotations)
snAnno_clean |> count(snType)

# using match to add annotated name columns for wt10
sce$snAnno <- snAnno_clean$snType[match(sce$wT_10_Erik, snAnno_clean$Cluster)]

# Now that sce has been established with new annotations from Excel sheet, we can
# plot against qc metrics (also in colData of sce with annotation identities):

####### Exploring Barcodes  #####################################################
UMI_ct_k <- ggcells(sce, mapping = aes(x = cellType_k, y = sum, fill = cellType_k)) +
  geom_boxplot() +
  scale_fill_manual(values = cell_type_colors[levels(sce$cellType_k)], drop = TRUE) +
  my_theme +
  theme(
    legend.position = "None", axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

