## May 17, 2023 - Bukola Ajanaku
# Updated UMAPs with confirmed color palette, finalizing axes titles, separating out 
# OPC_noisy class, and printing pre/post Harmony together!
# qrsh -l mem_free=40G,h_vmem=40G

library("SingleCellExperiment")
library("here")
library("tidyverse")
library("scater")
library("ggplot2")
library("cowplot")

# loading pre-harmony sce object 
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
          "sce_uncorrected_PCA.Rdata"), verbose = TRUE)
# sce_uncorrected

# Loading post-harmony sce object
load(file = here("processed-data", "99_paper_figs", "sce_objects", 
                 "sce_final_preHbdrop.RDATA"), verbose = TRUE)
# sce 

# sourcing official color palette 
source(file = here("code", "99_paper_figs", "source_colors.R"))
# bulk_colors and sn_colors

# creating plot directory 
plot_dir <- here("plots", "99_paper_figs", "02c_UMAPs_final")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

# keeping undropped sce as sce_dirty
sce_dirty <- sce 

# cleaning sce object
sce <- sce[, sce$final_Annotations != "OPC_noisy"]
sce <- sce[, sce$final_Annotations != "Excit.Neuron"]

# UMAP pre and post harmony by Sample 
pdf(here(plot_dir, "UMAP_PrePostHarmony_bySample.pdf"), 
    width = 11, height = 9)

plot1_unc <- plotReducedDim(sce_uncorrected, dimred = "UMAP") +
  geom_point(aes(color = sce_uncorrected$Sample)) + 
  theme(legend.position = "none") + 
  labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2")

plot2_unc <- plotReducedDim(sce_uncorrected, dimred = "UMAP") +
  geom_point(aes(color = sce_uncorrected$Sample)) + 
  facet_wrap(~ sce_uncorrected$Sample) + 
  guides(color = guide_legend(title="Sample ID")) + 
  labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2")

plot1 <- plotReducedDim(sce, dimred = "UMAP") +
  geom_point(aes(color = sce$Sample)) + 
  theme(legend.position = "none") + 
  labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2")

plot2 <- plotReducedDim(sce, dimred = "UMAP") +
  geom_point(aes(color = sce$Sample)) + 
  facet_wrap(~ sce$Sample) + 
  guides(color = guide_legend(title="Sample ID")) + 
  labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2")

plot_grid(plot1_unc, plot2_unc, plot1, plot2,
          nrow = 2,
          labels = c("A", "", "B", ""))
dev.off()

# UMAP pre and post harmony colored by Sample, faceted by Run
  # changing Run Labels 
sce_uncorrected$Run <- recode(sce_uncorrected$Run,
                              "1" = "Run 1",
                              "2" = "Run 2",
                              "3" = "Run 3")
sce$Run <- recode(sce$Run,
                  "1" = "Run 1",
                  "2" = "Run 2",
                  "3" = "Run 3")

pdf(here(plot_dir, "UMAP_PrePostHarmony_bySample_facetbyRun.pdf"), 
    width = 16, height = 9)

plot1_unc <- plotReducedDim(sce_uncorrected, dimred = "UMAP") +
  geom_point(aes(color = sce_uncorrected$Sample)) + 
  theme(legend.position = "none") + 
  labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2")

plot2_unc <- plotReducedDim(sce_uncorrected, dimred = "UMAP") +
  geom_point(aes(color = sce_uncorrected$Sample)) + 
  facet_wrap(~ sce_uncorrected$Run) + 
  guides(color = guide_legend(title="Sample ID")) + 
  labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2")

plot1 <- plotReducedDim(sce, dimred = "UMAP") +
  geom_point(aes(color = sce$Sample)) + 
  theme(legend.position = "none") + 
  labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2")

plot2 <- plotReducedDim(sce, dimred = "UMAP") +
  geom_point(aes(color = sce$Sample)) + 
  facet_wrap(~ sce$Run) + 
  guides(color = guide_legend(title="Sample ID")) + 
  labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2")

plot_grid(plot1_unc, plot2_unc, plot1, plot2,
          nrow = 2,
          labels = c("A", "", "B", ""),
          rel_widths = c(.3, .6))
dev.off()

## UMAP by cell type facet wrapped
pdf(here(plot_dir, "UMAP_harmony_by_CellType_faceted.pdf"),
    width = 15, height = 8)
plot1 <- plotReducedDim(sce, dimred = "UMAP") +
  geom_point(aes(color = sce$final_Annotations)) +
  scale_colour_manual(values = sn_colors) +
  theme(legend.position = "none") + 
  labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2")

plot2 <- plotReducedDim(sce, dimred = "UMAP") +
  geom_point(aes(color = sce$final_Annotations)) +
  scale_colour_manual(values = sn_colors) +
  facet_wrap(~ sce$final_Annotations) +
  guides(color = guide_legend(title="Cell Type")) + 
  labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2")

plot_grid(plot1, plot2)
dev.off()

##### for OneDrive Pieces ######################################################
# UMAP pre and post harmony by Sample 
# pre-harmony
pdf(here(plot_dir, "forOneDrive" , "sfigu_preHarmony_UMAP_bySample.pdf"), 
    width = 7, height = 3.5)

  plot1_unc <- plotReducedDim(sce_uncorrected, dimred = "UMAP") +
    geom_point(aes(color = sce_uncorrected$Sample)) + 
    theme(legend.position = "none") + 
    labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2")
  plot2_unc <- plotReducedDim(sce_uncorrected, dimred = "UMAP") +
    geom_point(aes(color = sce_uncorrected$Sample)) + 
    facet_wrap(~ sce_uncorrected$Sample) + 
    guides(color = guide_legend(title="Sample ID")) + 
    labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2")
      
      plot_grid(plot1_unc, plot2_unc,
                nrow = 1)
dev.off()

# post-harmony
pdf(here(plot_dir, "forOneDrive" , "sfigu_preHarmony_UMAP_bySample.pdf"), 
    width = 7, height = 3.5)

  plot1 <- plotReducedDim(sce, dimred = "UMAP") +
    geom_point(aes(color = sce$Sample)) + 
    theme(legend.position = "none") + 
    labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2")
  plot2 <- plotReducedDim(sce, dimred = "UMAP") +
    geom_point(aes(color = sce$Sample)) + 
    facet_wrap(~ sce$Sample) + 
    guides(color = guide_legend(title="Sample ID")) + 
    labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2")
  
      plot_grid(plot1, plot2,
              nrow = 1)
dev.off()

# UMAP pre and post harmony colored by Sample, faceted by Run
# changing Run Labels 
sce_uncorrected$Run <- recode(sce_uncorrected$Run,
                              "1" = "Run 1",
                              "2" = "Run 2",
                              "3" = "Run 3")
sce$Run <- recode(sce$Run,
                  "1" = "Run 1",
                  "2" = "Run 2",
                  "3" = "Run 3")

# pre-harmony
pdf(here(plot_dir, "forOneDrive" , "sfigu_preHarmony_UMAP_byRun.pdf"), 
    width = 8, height = 3)

  plot1_unc <- plotReducedDim(sce_uncorrected, dimred = "UMAP") +
    geom_point(aes(color = sce_uncorrected$Sample)) + 
    theme(legend.position = "none") + 
    labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2")
  plot2_unc <- plotReducedDim(sce_uncorrected, dimred = "UMAP") +
    geom_point(aes(color = sce_uncorrected$Sample)) + 
    facet_wrap(~ sce_uncorrected$Run) + 
    guides(color = guide_legend(title="Sample ID")) + 
    labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2")

    plot_grid(plot1_unc, plot2_unc,
              nrow = 1,
              rel_widths = c(.3, .7))
dev.off()

# post-harmony
pdf(here(plot_dir, "forOneDrive" , "sfigu_postHarmony_UMAP_byRun.pdf"), 
    width = 8, height = 3)
  
  plot1 <- plotReducedDim(sce, dimred = "UMAP") +
    geom_point(aes(color = sce$Sample)) + 
    theme(legend.position = "none") + 
    labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2")
  
  plot2 <- plotReducedDim(sce, dimred = "UMAP") +
    geom_point(aes(color = sce$Sample)) + 
    facet_wrap(~ sce$Run) + 
    guides(color = guide_legend(title="Sample ID")) + 
    labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2")
    
    plot_grid(plot1_unc, plot2_unc, 
              nrow = 1,
              rel_widths = c(.3, .7))
dev.off()

## UMAP by cell type facet wrapped
pdf(here(plot_dir, "forOneDrive", "sfigu_UMAP_byCellType_faceted.pdf"),
    width = 9, height = 4)

    plot1 <- plotReducedDim(sce, dimred = "UMAP") +
      geom_point(aes(color = sce$final_Annotations)) +
      scale_colour_manual(values = sn_colors) +
      theme(legend.position = "none") + 
      labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2")
    
    plot2 <- plotReducedDim(sce, dimred = "UMAP") +
      geom_point(aes(color = sce$final_Annotations)) +
      scale_colour_manual(values = sn_colors) +
      facet_wrap(~ sce$final_Annotations) +
      guides(color = guide_legend(title="Cell Type")) + 
      labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2")
    
  plot_grid(plot1, plot2,
            rel_widths = c(.4, .6))
dev.off()



# 