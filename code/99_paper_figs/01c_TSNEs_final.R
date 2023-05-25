## May 16, 2023 - Bukola Ajanaku
# Updated TSNEs with confirmed color palette, finalizing axes titles, separating out 
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

sce_dirty <- sce 

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


# TSNE pre and post harmony by Sample 
pdf(here(plot_dir, "TSNE_PrePostHarmony_bySample.pdf"), 
    width = 11, height = 9)

plot1_unc <- plotReducedDim(sce_uncorrected, dimred = "TSNE") +
  geom_point(aes(color = sce_uncorrected$Sample)) + 
  theme(legend.position = "none") + 
  labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")

plot2_unc <- plotReducedDim(sce_uncorrected, dimred = "TSNE") +
  geom_point(aes(color = sce_uncorrected$Sample)) + 
  facet_wrap(~ sce_uncorrected$Sample) + 
  guides(color = guide_legend(title="Sample ID")) + 
  labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")

plot1 <- plotReducedDim(sce, dimred = "TSNE") +
  geom_point(aes(color = sce$Sample)) + 
  theme(legend.position = "none") + 
  labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")

plot2 <- plotReducedDim(sce, dimred = "TSNE") +
  geom_point(aes(color = sce$Sample)) + 
  facet_wrap(~ sce$Sample) + 
  guides(color = guide_legend(title="Sample ID")) + 
  labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")

plot_grid(plot1_unc, plot2_unc, plot1, plot2,
          nrow = 2,
          labels = c("A", "", "B", ""))
dev.off()

# TSNE pre and post harmony colored by Sample, faceted by Run

  # changing Run Labels 
sce_uncorrected$Run <- recode(sce_uncorrected$Run,
                              "1" = "Run 1",
                              "2" = "Run 2",
                              "3" = "Run 3")
sce$Run <- recode(sce$Run,
                  "1" = "Run 1",
                  "2" = "Run 2",
                  "3" = "Run 3")

pdf(here(plot_dir, "TSNE_PrePostHarmony_bySample_facetbyRun.pdf"), 
    width = 16, height = 9)

plot1_unc <- plotReducedDim(sce_uncorrected, dimred = "TSNE") +
  geom_point(aes(color = sce_uncorrected$Sample)) + 
  theme(legend.position = "none") + 
  labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")

plot2_unc <- plotReducedDim(sce_uncorrected, dimred = "TSNE") +
  geom_point(aes(color = sce_uncorrected$Sample)) + 
  facet_wrap(~ sce_uncorrected$Run) + 
  guides(color = guide_legend(title="Sample ID")) + 
  labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")

plot1 <- plotReducedDim(sce, dimred = "TSNE") +
  geom_point(aes(color = sce$Sample)) + 
  theme(legend.position = "none") + 
  labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")

plot2 <- plotReducedDim(sce, dimred = "TSNE") +
  geom_point(aes(color = sce$Sample)) + 
  facet_wrap(~ sce$Run) + 
  guides(color = guide_legend(title="Sample ID")) + 
  labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")

plot_grid(plot1_unc, plot2_unc, plot1, plot2,
          nrow = 2,
          labels = c("A", "", "B", ""),
          rel_widths = c(.3, .6))
dev.off()

## TSNE by cell type
pdf(here(plot_dir, "TSNE_harmony_by_CellType.pdf"))
plotReducedDim(sce, dimred = "TSNE") +
  geom_point(aes(color = sce$final_Annotations)) + 
  scale_colour_manual(values = sn_colors) + 
  guides(color = guide_legend(title="Cell Type")) + 
  labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")
dev.off()

## TSNE by cell type facet wrapped
pdf(here(plot_dir, "TSNE_harmony_by_CellType_faceted.pdf"),
    width = 15, height = 8)
plot1 <- plotReducedDim(sce, dimred = "TSNE") +
  geom_point(aes(color = sce$final_Annotations)) +
  scale_colour_manual(values = sn_colors) +
  theme(legend.position = "none") + 
  labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")

plot2 <- plotReducedDim(sce, dimred = "TSNE") +
  geom_point(aes(color = sce$final_Annotations)) +
  scale_colour_manual(values = sn_colors) +
  facet_wrap(~ sce$final_Annotations) +
  guides(color = guide_legend(title="Cell Type")) + 
  labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")

plot_grid(plot1, plot2)
dev.off()

# TSNE by cell type pre and post drop
pdf(here(plot_dir, "TSNE_by_CellType_PrePostDrop.pdf"), 
    width = 11, height = 9)
  plot1_dirty <- plotReducedDim(sce_dirty, dimred = "TSNE") +
    geom_point(aes(color = sce_dirty$final_Annotations)) + 
    theme(legend.position = "none") +
    guides(color = guide_legend(title="Cell Type")) + 
    scale_colour_manual(values = sn_colors) + 
    labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")
  
  
  plot2_dirty <- plotReducedDim(sce_dirty, dimred = "TSNE") +
    geom_point(aes(color = sce_dirty$final_Annotations)) + 
    facet_wrap(~ sce_dirty$final_Annotations) + 
    guides(color = guide_legend(title="Cell Type")) + 
    scale_colour_manual(values = sn_colors) + 
    labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")
  
  
  plot1 <- plotReducedDim(sce, dimred = "TSNE") +
    geom_point(aes(color = sce$final_Annotations)) + 
    theme(legend.position = "none") +
    guides(color = guide_legend(title="Cell Type")) + 
    scale_colour_manual(values = sn_colors) + 
    labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")
  
  
  plot2 <- plotReducedDim(sce, dimred = "TSNE") +
    geom_point(aes(color = sce$final_Annotations)) + 
    facet_wrap(~ sce$final_Annotations) + 
    guides(color = guide_legend(title="Cell Type")) + 
    scale_colour_manual(values = sn_colors) + 
    labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")
  
  
  plot_grid(plot1_dirty, plot2_dirty, plot1, plot2,
            nrow = 2,
            labels = c("A", "", "B", ""))
dev.off()

# TSNE pre and post Harmony, colored by NeuN-sorting and faceted by Sample ID.
pdf(here(plot_dir, "TSNE_PrePostHarmony_by_NeuN.pdf"), width = 10, height = 9)

# grabbing NeuN data for sce uncorrected
  sce_uncorrected$NeuN <- sce_dirty$NeuN[match(sce_uncorrected$Barcode, sce_dirty$Barcode)]
  
      plot1_unc <- plotReducedDim(sce_uncorrected, dimred = "TSNE") +
        geom_point(aes(color = sce_uncorrected$NeuN)) + 
        theme(legend.position = "none") + 
        labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")
      
      plot2_unc <- plotReducedDim(sce_uncorrected, dimred = "TSNE") +
        geom_point(aes(color = sce_uncorrected$NeuN)) +
        facet_wrap(~ sce_uncorrected$Sample) + 
        scale_color_discrete(name = "NeuN Sorted?", labels=c('Yes', 'No')) + 
        labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")
      
      plot1 <- plotReducedDim(sce, dimred = "TSNE") +
        geom_point(aes(color = sce$NeuN)) + 
        theme(legend.position = "none") + 
        labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")
      
      plot2 <- plotReducedDim(sce, dimred = "TSNE") +
        geom_point(aes(color = sce$NeuN)) +
        facet_wrap(~ sce$Sample) + 
        scale_color_discrete(name = "NeuN Sorted?", labels=c('Yes', 'No')) + 
        labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")
      
      plot_grid(plot1_unc, plot2_unc, plot1, plot2,
                nrow = 2,
                labels = c("A", "", "B", ""))
dev.off()


###### FOR ONE DRIVE PIECES ####################################################
# 1 & 2:  TSNE pre and post harmony by Sample 

# pre-harmony
pdf(here(plot_dir, "forOneDrive","sfigu_preHarmony_TSNE_bySample.pdf"), 
    width = 7, height = 3.5)

  plot1_unc <- plotReducedDim(sce_uncorrected, dimred = "TSNE") +
    geom_point(aes(color = sce_uncorrected$Sample)) + 
    theme(legend.position = "none") + 
    labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")
  plot2_unc <- plotReducedDim(sce_uncorrected, dimred = "TSNE") +
    geom_point(aes(color = sce_uncorrected$Sample)) + 
    facet_wrap(~ sce_uncorrected$Sample) + 
    guides(color = guide_legend(title="Sample ID")) + 
    labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")


    plot_grid(plot1_unc, plot2_unc,
              nrow = 1)
dev.off()

# post-harmoy 
pdf(here(plot_dir, "forOneDrive", "sfigu_postHarmony_TSNE_bySample.pdf"), 
    width = 7, height = 3.5)

  plot1 <- plotReducedDim(sce, dimred = "TSNE") +
    geom_point(aes(color = sce$Sample)) + 
    theme(legend.position = "none") + 
    labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")
  plot2 <- plotReducedDim(sce, dimred = "TSNE") +
    geom_point(aes(color = sce$Sample)) + 
    facet_wrap(~ sce$Sample) + 
    guides(color = guide_legend(title="Sample ID")) + 
    labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")

    plot_grid(plot1, plot2,
              nrow = 1)
dev.off()


# 3 & 4: TSNE pre and post harmony colored by Sample, faceted by Run
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
pdf(here(plot_dir, "forOneDrive", "sfigu_preHarmony_TSNE_byRun.pdf"), 
    width = 8, height = 3)

  plot1_unc <- plotReducedDim(sce_uncorrected, dimred = "TSNE") +
    geom_point(aes(color = sce_uncorrected$Sample)) + 
    theme(legend.position = "none") + 
    labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")
  plot2_unc <- plotReducedDim(sce_uncorrected, dimred = "TSNE") +
    geom_point(aes(color = sce_uncorrected$Sample)) + 
    facet_wrap(~ sce_uncorrected$Run) + 
    guides(color = guide_legend(title="Sample ID")) + 
    labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")

    plot_grid(plot1_unc, plot2_unc,
              nrow = 1,
              rel_widths = c(.3, .7))
dev.off()

# post-harmony 
pdf(here(plot_dir, "forOneDrive", "sfigu_postHarmony_TSNE_byRun.pdf"), 
    width = 8, height = 3)

  plot1 <- plotReducedDim(sce, dimred = "TSNE") +
    geom_point(aes(color = sce$Sample)) + 
    theme(legend.position = "none") + 
    labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")
  plot2 <- plotReducedDim(sce, dimred = "TSNE") +
    geom_point(aes(color = sce$Sample)) + 
    facet_wrap(~ sce$Run) + 
    guides(color = guide_legend(title="Sample ID")) + 
    labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")

    plot_grid(plot1, plot2,
              nrow = 1,
              rel_widths = c(.3, .7))
dev.off()


## TSNE by cell type
pdf(here(plot_dir, "forOneDrive", "mfigu_TSNE_byCellType_noFacet.pdf"))
  plotReducedDim(sce, dimred = "TSNE") +
    geom_point(aes(color = sce$final_Annotations)) + 
    scale_colour_manual(values = sn_colors) + 
    guides(color = guide_legend(title="Cell Type")) + 
    labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")
dev.off()

## TSNE by cell type facet wrapped
pdf(here(plot_dir,  "forOneDrive", "mfigu_TSNE_byCellType_faceted.pdf"),
    width = 7, height = 4)
    plot1 <- plotReducedDim(sce, dimred = "TSNE") +
      geom_point(aes(color = sce$final_Annotations)) +
      scale_colour_manual(values = sn_colors) +
      theme(legend.position = "none") + 
      labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")
    plot2 <- plotReducedDim(sce, dimred = "TSNE") +
      geom_point(aes(color = sce$final_Annotations)) +
      scale_colour_manual(values = sn_colors) +
      facet_wrap(~ sce$final_Annotations) +
      guides(color = guide_legend(title="Cell Type")) + 
      labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")

      plot_grid(plot1, plot2)
dev.off()


# TSNE by cell type pre and post drop (already harmonized)
# pre-drop 
pdf(here(plot_dir, "forOneDrive", "sfigu_TSNE_harmony_preDrop_byCellType.pdf"), 
    width = 7, height = 4)
    plot1_dirty <- plotReducedDim(sce_dirty, dimred = "TSNE") +
      geom_point(aes(color = sce_dirty$final_Annotations)) + 
      theme(legend.position = "none") +
      guides(color = guide_legend(title="Cell Type")) + 
      scale_colour_manual(values = sn_colors) + 
      labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")
    
    plot2_dirty <- plotReducedDim(sce_dirty, dimred = "TSNE") +
      geom_point(aes(color = sce_dirty$final_Annotations)) + 
      facet_wrap(~ sce_dirty$final_Annotations) + 
      guides(color = guide_legend(title="Cell Type")) + 
      scale_colour_manual(values = sn_colors) + 
      labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")

      
      plot_grid(plot1_dirty, plot2_dirty,
                nrow = 1)
dev.off()

# post-drop
pdf(here(plot_dir, "forOneDrive", "sfigu_TSNE_harmony_postDrop_byCellType.pdf"), 
    width = 7, height = 4)

  plot1 <- plotReducedDim(sce, dimred = "TSNE") +
    geom_point(aes(color = sce$final_Annotations)) + 
    theme(legend.position = "none") +
    guides(color = guide_legend(title="Cell Type")) + 
    scale_colour_manual(values = sn_colors) + 
    labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")
  
  
  plot2 <- plotReducedDim(sce, dimred = "TSNE") +
    geom_point(aes(color = sce$final_Annotations)) + 
    facet_wrap(~ sce$final_Annotations) + 
    guides(color = guide_legend(title="Cell Type")) + 
    scale_colour_manual(values = sn_colors) + 
    labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")
  
    plot_grid(plot1, plot2,
              nrow = 1)
dev.off()


# TSNE pre and post Harmony, colored by NeuN-sorting and faceted by Sample ID.
# grabbing NeuN data for sce uncorrected
sce_uncorrected$NeuN <- sce_dirty$NeuN[match(sce_uncorrected$Barcode, sce_dirty$Barcode)]

# pre-harmony
pdf(here(plot_dir, "forOneDrive", "sfigu_TSNE_preHarmony_colorByNeuNSort.pdf"), 
    width = 7, height = 3.5)

  plot1_unc <- plotReducedDim(sce_uncorrected, dimred = "TSNE") +
    geom_point(aes(color = sce_uncorrected$NeuN)) + 
    theme(legend.position = "none") + 
    labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")
  plot2_unc <- plotReducedDim(sce_uncorrected, dimred = "TSNE") +
    geom_point(aes(color = sce_uncorrected$NeuN)) +
    facet_wrap(~ sce_uncorrected$Sample) + 
    scale_color_discrete(name = "NeuN Sorted?", labels=c('Yes', 'No')) + 
    labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")
    
    plot_grid(plot1_unc, plot2_unc, 
              nrow = 1)
dev.off()

# post-harmony 
pdf(here(plot_dir, "forOneDrive", "sfigu_TSNE_postHarmony_colorByNeuNSort.pdf"), 
    width = 7, height = 3.5)
  plot1 <- plotReducedDim(sce, dimred = "TSNE") +
    geom_point(aes(color = sce$NeuN)) + 
    theme(legend.position = "none") + 
    labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")
  
  plot2 <- plotReducedDim(sce, dimred = "TSNE") +
    geom_point(aes(color = sce$NeuN)) +
    facet_wrap(~ sce$Sample) + 
    scale_color_discrete(name = "NeuN Sorted?", labels=c('Yes', 'No')) + 
    labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")

    plot_grid(plot1, plot2,
              nrow = 1)
dev.off()

# Done.
