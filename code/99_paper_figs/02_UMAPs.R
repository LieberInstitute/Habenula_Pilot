## April 6, 2023 - Bukola Ajanaku
# UMAPs by my celltypes and custom color pallete 
# qrsh -l mem_free=20G,h_vmem=20G

library("SingleCellExperiment")
library("here")
library("tidyverse")
library("scater")
library("ggplot2")
library("cowplot")


# Loading sce object
load(here("processed-data", "99_paper_figs",  "sce_objects", "sce_final_preHbdrop.RDATA"),
     verbose = TRUE)
  # sce_final_preHbdrop , sce_sorted, sce_unsorted

table(sce_final_preHbdrop$final_Annotations)
  # Astrocyte         Endo Excit.Neuron   Excit.Thal   Inhib.Thal        LHb.1 
  # 538           38           51         1800         7612          201 
  # LHb.2        LHb.3        LHb.4        LHb.5        LHb.6        LHb.7 
  # 266          134          477           83           39         1014 
  # MHb.1        MHb.2        MHb.3    Microglia        Oligo          OPC 
  # 152          540           18          145         2178         1796 

# creating plot directory
plot_dir <- here("plots", "99_paper_figs", "02_UMAPs")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

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

##### PLOTTING UMAPs ###########################################################
pdf(here(plot_dir, "UMAP_harmony_by_finalAnno.pdf"))
  plotReducedDim(sce_final_preHbdrop, dimred = "UMAP") +
    geom_point(aes(color = sce_final_preHbdrop$final_Annotations), alpha = 0.4) + 
    scale_colour_manual(values = cluster_colors) + 
    guides(color = guide_legend(title="Cell Type"))
dev.off()

pdf(here(plot_dir, "UMAP_harmony_by_final_Annotations_splitbyfinalAnno.pdf"),
    width = 12, height = 7)
  plot1 <- plotReducedDim(sce_final_preHbdrop, dimred = "UMAP") +
    geom_point(aes(color = sce_final_preHbdrop$final_Annotations), alpha = 0.4) + 
    scale_colour_manual(values = cluster_colors) + 
    theme(legend.position = "none")

  plot2 <- plotReducedDim(sce_final_preHbdrop, dimred = "UMAP") +
    geom_point(aes(color = sce_final_preHbdrop$final_Annotations), alpha = 0.4) + 
    scale_colour_manual(values = cluster_colors) +
    facet_wrap(~ sce_final_preHbdrop$final_Annotations) + 
    guides(color = guide_legend(title="Cell Type"))
  
  plot_grid(plot1, plot2)
dev.off()

# 
plot_sorted <-   plotReducedDim(sce_sorted, dimred = "UMAP") +
  geom_point(aes(color = sce_sorted$final_Annotations), alpha = 0.4) + 
  scale_colour_manual(values = cluster_colors) +
  facet_grid(sce_sorted$NeuN ~ sce_sorted$Sample) + 
  guides(color = guide_legend(title="Cell Type"))

plot_unsorted <-   plotReducedDim(sce_unsorted, dimred = "UMAP") +
  geom_point(aes(color = sce_unsorted$final_Annotations), alpha = 0.4) + 
  scale_colour_manual(values = cluster_colors) +
  facet_grid(sce_unsorted$NeuN ~ sce_unsorted$Sample) + 
  guides(color = guide_legend(title="Cell Type"))


pdf(here(plot_dir, "UMAP_harmony_by_finalAnno_splitbySampleAndSorting.pdf"), width = 13, height = 9)
plot_grid(
  plot_sorted,
  plot_unsorted,
  ncol = 1
)
dev.off()


# 