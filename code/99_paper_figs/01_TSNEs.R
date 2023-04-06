## April 6, 2023 - Bukola Ajanaku
# TSNEs by my celltypes and custom color pallete 
# qrsh -l mem_free=20G,h_vmem=20G

library("SingleCellExperiment")
library("here")
library("tidyverse")
library("scater")
library("ggplot2")
library("cowplot")


# Loading sce object with snAnno2 (Has combined MHb clusters and adjusted LHb clusters. 
# We probs wanna rename the general "Hb" cluster to something like Exc.Neuron!!)
load(here("processed-data", "04_snRNA-seq", "sce_objects", "sce_with_snAnno2.RDATA"), 
     verbose = TRUE)

# creating plot directory
plot_dir <- here("plots", "99_paper_figs", "01_TSNEs")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

# checking sce object
table(colData(sce$snAnno2))
  # Astrocyte       Endo Excit.Thal  Hb Inhib.Thal      LHb.1      LHb.2 
  # 538         38       1800         51       7612        201        266 
  # LHb.3      LHb.4      LHb.5      LHb.7      LHb.8      MHb.1      MHb.2 
  # 134        477         83         39       1014        152        145 
  # MHb.3      MHb.4  Microglia      Oligo        OPC 
  # 395         18        145       2178       1796 

# changing general "Hb" cluster to Excit.Neuron
sce$snAnno2[sce$snAnno2 == "Hb"] <- "Excit.Neuron"

# adding the NeuN sorting/unsorting column to sce 
sce$NeuN <- "NeuN.Unsorted"
sce$NeuN[sce$Sample %in% c("Br1092", "Br1204", "Br5555", "Br5558")] <- "NeuN.Sorted"

# checking
table(sce$Sample, sce$NeuN)
#         NeuN.Sorted NeuN.Unsorted
# Br1092        3212             0
# Br1204        1239             0
# Br1469           0          2151
# Br1735           0          3101
# Br5555        3542             0
# Br5558         778             0
# Br5639           0          3059


# adding color pallete (same color scheme used for progress report heatmap)
cluster_colors <- c( "Oligo" = c("#5C4033"), # dark grown
                       "OPC"= c("#899499"), # pewter (grey)
                       "Microglia" = c("#4CBB17"), # kelly green
                       "Astrocyte" = c("#CC5500"), # burnt orange
                       "Endo" = c("#702963"), # byzantium
                       "Excit.Neuron" = c("#FAFA33"), # lemon yellow
                       "Inhib.Thal" = c('#FF0000'), #red
                       "Excit.Thal" = c('#0818A8'), # zaffre (dark blue)
                       "LHb.1" = c("#5F9EA0"), # cadet Blue
                       "LHb.2" = c("#5D3FD3"), # iris
                       "LHb.3" = c ("#4682B4"), #Steel Blue
                       "LHb.4" = c("#1F51FF"), # neon blue
                       "LHb.5" = c("#6495ED"), # Cornflower Blue
                       "LHb.6" = c("#088F8F"), # Blue Green 
                       "LHb.7" = c("#4169E1"), # royal bluee
                       "MHb.1" = c("#40E0D0"), # Turquoise
                       "MHb.2" = c("#96DED1"), #Robin Egg Blue
                       "MHb.3" = c("#7DF9FF") # light blue
)

##### PLOTTING TSNEs ###########################################################
pdf(here(plot_dir, "TSNE_harmony_by_snAnno2.pdf"))
  plotReducedDim(sce, dimred = "TSNE", colour_by = "snAnno2") +
    scale_colour_manual(values = cluster_colors)
dev.off()

pdf(here(plot_dir, "TSNE_harmony_by_snAnno2_splitbysnAnno2.pdf"))
  plotReducedDim(sce, dimred = "TSNE", colour_by = "snAnno2") +
    scale_colour_manual(values = cluster_colors) +
    facet_wrap(~ sce$snAnno2)
dev.off()

# Separating sce by NeuN sorting
sce_sorted <- sce[, sce$NeuN == "NeuN.Sorted"] 
sce_unsorted <- sce[, sce$NeuN == "NeuN.Unsorted"] 

# 
plot_sorted <- plotReducedDim(sce_sorted, dimred = "TSNE", colour_by = "snAnno2") +
  scale_colour_manual(values = cluster_colors) +
  facet_grid(sce_sorted$NeuN ~ sce_sorted$Sample)

plot_unsorted <- plotReducedDim(sce_unsorted, dimred = "TSNE", colour_by = "snAnno2") +
  scale_colour_manual(values = cluster_colors) +
  facet_grid(sce_unsorted$NeuN ~ sce_unsorted$Sample)


pdf(here(plot_dir, "TSNE_harmony_by_snAnno2_splitbySampleAndSorting.pdf"), width = 10, height = 9)
  plot_grid(
    plot_sorted,
    plot_unsorted,
    ncol = 1
  )
dev.off()


# 
