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
library("jaffelab")

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
         Cluster = paste0("10wTrap_", Cluster)) |>
  mutate(splitProbClusts = paste(problemClusts, ss(Cluster, "_", 2), sep = "_")) |>
  mutate(splitSNType = paste(snType, ss(Cluster, "_", 2), sep = "_"))

# dropping three extra columns 
snAnno_clean <- select(snAnno_clean, -c(starts_with("Na"), "MoreInfo"))

# sanity check (checking for repeats in annotations to check that split by cluster makes sense)
snAnno_clean |> count(snType)
snAnno_clean |> count(problemClusts)

# using match to add annotated name columns for wt10
## actual annotatios for now  
  sce$snAnno <- snAnno_clean$snType[match(sce$wT_10_Erik, snAnno_clean$Cluster)]

## adding columns that we want to plot (annotations but separated by clusters and 
## annotations with problematic cluster highlights all split by clusters)
  sce$splitProbClusts <- snAnno_clean$splitProbClusts[match(sce$wT_10_Erik, snAnno_clean$Cluster)]
  sce$splitSNType <- snAnno_clean$splitSNType[match(sce$wT_10_Erik, snAnno_clean$Cluster)]
  
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
pdf(file = here(plot_dir, "wt10_LIBSIZE_QC.pdf"), height = 7, width = 11)
  # coloring by annotations with problems highlighted
  ggcells(sce, mapping = aes(x = splitProbClusts, y = sum, fill = splitProbClusts)) +
      geom_boxplot() +
      scale_fill_manual(values = grabColors(length(unique(sce$splitProbClusts)))) +
      theme_linedraw() +
      theme(
        legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
      ggtitle("Library Size with Problematic Cluster Highlights")
  
  # coloring by annotations only (not split)
  ggcells(sce, mapping = aes(x = snAnno, y = sum, fill = snAnno)) +
    geom_boxplot() +
    scale_fill_manual(values = grabColors(length(unique(sce$snAnno)))) +
    theme_linedraw() +
    theme(
      legend.position = "none", axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
    ggtitle("Library Size for all Annotated Clusters")
  
  # regular annotations split by cluster 
  ggcells(sce, mapping = aes(x = splitSNType, y = sum, fill = splitSNType)) +
    geom_boxplot() +
    scale_fill_manual(values = grabColors(length(unique(sce$splitSNType)))) +
    theme_linedraw() +
    theme(
      legend.position = "none", axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle("Library Size for Regular Annotations Split by Cluster")
dev.off()

# Checking detected features ("detected": the higher the better, dropping lows)
pdf(file = here(plot_dir, "wt10_DETECTED_QC.pdf"), height = 7, width = 11)
  # coloring by annotations with problems highlighted
  ggcells(sce, mapping = aes(x = splitProbClusts, y = detected, fill = splitProbClusts)) +
    geom_boxplot() +
    scale_fill_manual(values = grabColors(length(unique(sce$splitProbClusts)))) +
    theme_linedraw() +
    theme(
      legend.position = "none", axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle("Detected Features with Problematic Cluster Highlights")
  
  # coloring by annotations only (not split)
  ggcells(sce, mapping = aes(x = snAnno, y = detected, fill = snAnno)) +
    geom_boxplot() +
    scale_fill_manual(values = grabColors(length(unique(sce$snAnno)))) +
    theme_linedraw() +
    theme(
      legend.position = "none", axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle("Detected Features for all Annotated Clusters")
  
  # regular annotations split by cluster 
  ggcells(sce, mapping = aes(x = splitSNType, y = detected, fill = splitSNType)) +
    geom_boxplot() +
    scale_fill_manual(values = grabColors(length(unique(sce$splitSNType)))) +
    theme_linedraw() +
    theme(
      legend.position = "none", axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle("Detected Features for Regular Annotations Split by Cluster")
  
dev.off()

# Checking mitochondrial rate ("subsets_Mito_percent": the lower the better, dropping highs)
pdf(file = here(plot_dir, "wt10_MITO_PERCENT_QC.pdf"), height = 7, width = 11)
# coloring by annotations with problems highlighted
    ggcells(sce, mapping = aes(x = splitProbClusts, y = subsets_Mito_percent, fill = splitProbClusts)) +
      geom_boxplot() +
      scale_fill_manual(values = grabColors(length(unique(sce$splitProbClusts)))) +
      theme_linedraw() +
      theme(
        legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      ggtitle("Mito Percent with Problematic Cluster Highlights")
    
    # coloring by annotations only (not split)
    ggcells(sce, mapping = aes(x = snAnno, y = subsets_Mito_percent, fill = snAnno)) +
      geom_boxplot() +
      scale_fill_manual(values = grabColors(length(unique(sce$snAnno)))) +
      theme_linedraw() +
      theme(
        legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      ggtitle("Mito Percent for all Annotated Clusters")
    
    # regular annotations split by cluster 
    ggcells(sce, mapping = aes(x = splitSNType, y = subsets_Mito_percent, fill = splitSNType)) +
      geom_boxplot() +
      scale_fill_manual(values = grabColors(length(unique(sce$splitSNType)))) +
      theme_linedraw() +
      theme(
        legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      ggtitle("Mito Percent for Regular Annotations Split by Cluster")
dev.off()

# Checking doublelt scores ("doubletScore": the lower the better, typically dropping anything over 5)
pdf(file = here(plot_dir, "wt10_DOUBLET_SCORE_QC.pdf"), height = 7, width = 11)
# coloring by annotations with problems highlighted
    ggcells(sce, mapping = aes(x = splitProbClusts, y = doubletScore, fill = splitProbClusts)) +
      geom_boxplot() +
      scale_fill_manual(values = grabColors(length(unique(sce$splitProbClusts)))) +
      theme_linedraw() +
      theme(
        legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      ggtitle("Doublet Score with Problematic Cluster Highlights")
    
    # coloring by annotations only (not split)
    ggcells(sce, mapping = aes(x = snAnno, y = doubletScore, fill = snAnno)) +
      geom_boxplot() +
      scale_fill_manual(values = grabColors(length(unique(sce$snAnno)))) +
      theme_linedraw() +
      theme(
        legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      ggtitle("Doublet Score for all Annotated Clusters")
    
    # regular annotations split by cluster 
    ggcells(sce, mapping = aes(x = splitSNType, y = doubletScore, fill = splitSNType)) +
      geom_boxplot() +
      scale_fill_manual(values = grabColors(length(unique(sce$splitSNType)))) +
      theme_linedraw() +
      theme(
        legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      ggtitle("Doublet Score for Regular Annotations Split by Cluster")
dev.off()


# saving sce object
save(sce, file = here("processed-data", "04_snRNA-seq", "sce_objects", 
                                "sce_post_09_clustered_qc.Rdata"))
