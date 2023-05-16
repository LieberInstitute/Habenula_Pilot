## May 16, 2023 - Bukola Ajanaku
# Updated TSNEs with confirmed color pallete, finalizing axes titles, separating out 
# OPC_noisy class, and printing pre/post Harmony together!
# qrsh -l mem_free=40G,h_vmem=40G

library("SingleCellExperiment")
library("here")
library("tidyverse")
library("scater")
library("ggplot2")
library("cowplot")


# Loading sce object
load(file = here("processed-data", "99_paper_figs", "sce_objects", 
                 "sce_final_preHbdrop.RDATA"), verbose = TRUE)

# sourcing official color palette 
source(file = here("code", "99_paper_figs", "source_colors.R"))
  # bulk_colors and sn_colors

# creating plot directory 
plot_dir <- here("plots", "99_paper_figs", "01c_TSNEs_final")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

table(sce$final_Annotations)
  # Astrocyte         Endo Excit.Neuron   Excit.Thal   Inhib.Thal        LHb.1 
  # 538           38           51         1800         7612          201 
  # LHb.2        LHb.3        LHb.4        LHb.5        LHb.6        LHb.7 
  # 266          134          477           83           39         1014 
  # MHb.1        MHb.2        MHb.3    Microglia        Oligo          OPC 
  # 152          540           18          145         2178         1202 
  # OPC_noisy 
  # 594

sce_old <- sce 

# cleaning sce object
sce <- sce[, sce$final_Annotations != "OPC_noisy"]
sce <- sce[, sce$final_Annotations != "Excit.Neuron"]

table(sce$final_Annotations)
  # Astrocyte       Endo Excit.Thal Inhib.Thal      LHb.1      LHb.2      LHb.3 
  # 538         38       1800       7612        201        266        134 
  # LHb.4      LHb.5      LHb.6      LHb.7      MHb.1      MHb.2      MHb.3 
  # 477         83         39       1014        152        540         18 
  # Microglia      Oligo        OPC 
  # 145       2178       1202 


pdf(here(plot_dir, "TSNE_PrePostHarmony_bySample.pdf"), 
    width = 14, height = 9)

    plot1_old <- plotReducedDim(sce_old, dimred = "TSNE") +
      geom_point(aes(color = sce_old$Sample)) + 
      theme(legend.position = "none")
    
    plot2_old <- plotReducedDim(sce_old, dimred = "TSNE") +
      geom_point(aes(color = sce_old$Sample)) + 
      facet_wrap(~ sce_old$Sample) + 
      guides(color = guide_legend(title="Sample ID"))
    
    plot1 <- plotReducedDim(sce, dimred = "TSNE") +
      geom_point(aes(color = sce$Sample)) + 
      theme(legend.position = "none")
    
    plot2 <- plotReducedDim(sce, dimred = "TSNE") +
      geom_point(aes(color = sce$Sample)) + 
      facet_wrap(~ sce$Sample) + 
      guides(color = guide_legend(title="Sample ID"))

  plot_grid(plot1_old, plot2_old, plot1, plot2,
            nrow = 2)
dev.off()










# 