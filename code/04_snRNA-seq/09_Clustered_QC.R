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
library("scater")


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

# sourcing for custom  color palette 
source(here("code", "04_snRNA-seq", "sourcing", "color_Scheme_CT.R"))

# Now that sce has been established with new annotations from Excel sheet, we can
# plot against qc metrics (also in colData of sce with annotation identities):
####### Exploring Barcodes  #####################################################
# always create plot dir before plotting 
  plot_dir <- here("plots", "04_snRNA-seq", "09_Clustered_QC")
  if(!dir.exists(plot_dir)){
    dir.create(plot_dir)
  } 


# Checking library size ("sum": the higher the better, dropping lows)
pdf(file = here(plot_dir, "wt10_LIBSIZE_QC.pdf") )
  # coloring by walktrap cluster number
  ggcells(sce, mapping = aes(x = wT_10_Erik, y = sum, fill = wT_10_Erik)) +
    geom_boxplot() +
    scale_fill_manual(values = grabColors(length(unique(sce$wT_10_Erik)))) +
    theme_linedraw() +
    theme(
      legend.position = "Right" , axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1)
  )
   
  # coloring by snAnnotation clusters  
  ggcells(sce, mapping = aes(x = snAnno, y = sum, fill = snAnno)) +
      geom_boxplot() +
      scale_fill_manual(values = grabColors(length(unique(sce$snAnno)))) +
      theme_linedraw() +
      theme(
        legend.position = "Right", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)
  )
dev.off()


# Checking detected features ("detected": the higher the better, dropping lows)

# Checking mitochondrial rate ("subsets_Mito_percent": the lower the better, dropping highs)

# Checking doublelt scores ("doubletScore": the lower the better, typically dropping anything over 5)




