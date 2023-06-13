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
    # pdf
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
    # png
    png(here(plot_dir, "forOneDrive","sfigu_preHarmony_TSNE_bySample.png"), 
        width = 7, height = 3.5, units = "in", res = 1200)
    
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

# post-harmony 
    # pdf
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
    # png
    png(here(plot_dir, "forOneDrive", "sfigu_postHarmony_TSNE_bySample.png"), 
        width = 7, height = 3.5, units = "in", res = 1200)
    
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
# pdf
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
# png
png(here(plot_dir, "forOneDrive", "sfigu_preHarmony_TSNE_byRun.png"), 
    width = 8, height = 3, units = "in", res = 1200)

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

######

# post-harmony 
    # pdf
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
    # png
    png(here(plot_dir, "forOneDrive", "sfigu_postHarmony_TSNE_byRun.png"), 
        width = 8, height = 3, units = "in", res = 1200)
    
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

######

## TSNE by cell type
    # pdf
    pdf(here(plot_dir, "forOneDrive", "mfigu_TSNE_byCellType_noFacet.pdf"))
      plotReducedDim(sce, dimred = "TSNE") +
        geom_point(aes(color = sce$final_Annotations)) + 
        scale_colour_manual(values = sn_colors) + 
        guides(color = guide_legend(title="Cell Type")) + 
        labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")
    dev.off()
    # png
    png(here(plot_dir, "forOneDrive", "mfigu_TSNE_byCellType_noFacet.png"))
      plotReducedDim(sce, dimred = "TSNE") +
        geom_point(aes(color = sce$final_Annotations)) + 
        scale_colour_manual(values = sn_colors) + 
        guides(color = guide_legend(title="Cell Type")) + 
        labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")
    dev.off()

## TSNE by cell type facet wrapped
    # pdf 
    pdf(here(plot_dir,  "forOneDrive", "mfigu_TSNE_byCellType_faceted.pdf"),
        width = 9, height = 4)
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
    
          plot_grid(plot1, plot2,
                    rel_widths = c(.4, .6))
    dev.off()
    # png
    png(here(plot_dir,  "forOneDrive", "mfigu_TSNE_byCellType_faceted.png"),
        width = 9, height = 4, units = "in", res = 1200)
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
    
    plot_grid(plot1, plot2,
              rel_widths = c(.4, .6))
    dev.off()


# TSNE by cell type pre and post drop (already harmonized)
      # pre-drop 
        # pdf
      pdf(here(plot_dir, "forOneDrive", "sfigu_TSNE_harmony_preDrop_byCellType.pdf"), 
          width = 9, height = 4)
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
                      nrow = 1,
                      rel_widths = c(.4, .6))
      dev.off()
        # png
      png(here(plot_dir, "forOneDrive", "sfigu_TSNE_harmony_preDrop_byCellType.png"), 
          width = 9, height = 4, units = "in", res = 1200)
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
                nrow = 1,
                rel_widths = c(.4, .6))
      dev.off()

  # post-drop
    # pdf 
    pdf(here(plot_dir, "forOneDrive", "sfigu_TSNE_harmony_postDrop_byCellType.pdf"), 
        width = 9, height = 4)
    
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
                  nrow = 1,
                  rel_widths = c(.4, .6))
    dev.off()
    # png
    png(here(plot_dir, "forOneDrive", "sfigu_TSNE_harmony_postDrop_byCellType.png"), 
        width = 9, height = 4, units = "in", res = 1200)
    
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
              nrow = 1,
              rel_widths = c(.4, .6))
    dev.off()

# TSNE pre and post Harmony, colored by NeuN-sorting and faceted by Sample ID.
# grabbing NeuN data for sce uncorrected
sce_uncorrected$NeuN <- sce_dirty$NeuN[match(sce_uncorrected$Barcode, sce_dirty$Barcode)]

# pre-harmony
    # pdf 
    pdf(here(plot_dir, "forOneDrive", "sfigu_TSNE_preHarmony_colorByNeuNSort.pdf"), 
        width = 7.5, height = 3.5)
    
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
                  nrow = 1,
                  rel_widths = c(.4, .6))
    dev.off()
    # png
    png(here(plot_dir, "forOneDrive", "sfigu_TSNE_preHarmony_colorByNeuNSort.png"), 
        width = 7.5, height = 3.5, units = "in", res = 1200)
    
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
              nrow = 1,
              rel_widths = c(.4, .6))
    dev.off()

    # post-harmony 
    # pdf 
    pdf(here(plot_dir, "forOneDrive", "sfigu_TSNE_postHarmony_colorByNeuNSort.pdf"), 
        width = 7.5, height = 3.5)
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
                  nrow = 1,
                  rel_widths = c(.4, .6))
    dev.off()
    # png 
    png(here(plot_dir, "forOneDrive", "sfigu_TSNE_postHarmony_colorByNeuNSort.png"), 
        width = 7.5, height = 3.5, units = "in", res = 1200)
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
              nrow = 1,
              rel_widths = c(.4, .6))
    dev.off()


sessioninfo::session_info()

# ─ Session info ──────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.2.3 Patched (2023-04-07 r84211)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2023-06-13
# pandoc   2.19.2 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/bin/pandoc
# 
# ─ Packages ──────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# beachmat               2.14.2    2023-04-07 [2] Bioconductor
# beeswarm               0.4.0     2021-06-01 [2] CRAN (R 4.2.1)
# Biobase              * 2.58.0    2022-11-01 [2] Bioconductor
# BiocGenerics         * 0.44.0    2022-11-01 [2] Bioconductor
# BiocNeighbors          1.16.0    2022-11-01 [2] Bioconductor
# BiocParallel           1.32.6    2023-03-17 [2] Bioconductor
# BiocSingular           1.14.0    2022-11-01 [2] Bioconductor
# bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.2.1)
# cli                    3.6.1     2023-03-23 [2] CRAN (R 4.2.3)
# codetools              0.2-19    2023-02-01 [3] CRAN (R 4.2.3)
# colorout             * 1.2-2     2023-02-09 [1] Github (jalvesaq/colorout@79931fd)
# colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.2.2)
# cowplot              * 1.1.1     2020-12-30 [2] CRAN (R 4.2.1)
# DelayedArray           0.24.0    2022-11-01 [2] Bioconductor
# DelayedMatrixStats     1.20.0    2022-11-01 [2] Bioconductor
# dplyr                * 1.1.1     2023-03-22 [2] CRAN (R 4.2.3)
# fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.2.2)
# forcats              * 1.0.0     2023-01-29 [2] CRAN (R 4.2.2)
# generics               0.1.3     2022-07-05 [2] CRAN (R 4.2.1)
# GenomeInfoDb         * 1.34.9    2023-02-02 [2] Bioconductor
# GenomeInfoDbData       1.2.9     2022-09-29 [2] Bioconductor
# GenomicRanges        * 1.50.2    2022-12-16 [2] Bioconductor
# ggbeeswarm             0.7.1     2022-12-16 [2] CRAN (R 4.2.2)
# ggplot2              * 3.4.2     2023-04-03 [2] CRAN (R 4.2.3)
# ggrepel                0.9.3     2023-02-03 [2] CRAN (R 4.2.2)
# glue                   1.6.2     2022-02-24 [2] CRAN (R 4.2.1)
# gridExtra              2.3       2017-09-09 [2] CRAN (R 4.2.1)
# gtable                 0.3.3     2023-03-21 [2] CRAN (R 4.2.3)
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.2.1)
# hms                    1.1.3     2023-03-21 [2] CRAN (R 4.2.3)
# IRanges              * 2.32.0    2022-11-01 [2] Bioconductor
# irlba                  2.3.5.1   2022-10-03 [2] CRAN (R 4.2.1)
# lattice                0.20-45   2021-09-22 [3] CRAN (R 4.2.3)
# lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.2.1)
# lubridate            * 1.9.2     2023-02-10 [2] CRAN (R 4.2.2)
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.2.1)
# Matrix                 1.5-4     2023-04-04 [2] CRAN (R 4.2.3)
# MatrixGenerics       * 1.10.0    2022-11-01 [2] Bioconductor
# matrixStats          * 1.0.0     2023-06-02 [1] CRAN (R 4.2.3)
# munsell                0.5.0     2018-06-12 [2] CRAN (R 4.2.1)
# pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.2.3)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.2.1)
# purrr                * 1.0.1     2023-01-10 [2] CRAN (R 4.2.2)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.2.1)
# Rcpp                   1.0.10    2023-01-22 [2] CRAN (R 4.2.2)
# RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.2.3)
# readr                * 2.1.4     2023-02-10 [2] CRAN (R 4.2.2)
# rlang                  1.1.1     2023-04-28 [1] CRAN (R 4.2.3)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.2.1)
# rsvd                   1.0.5     2021-04-16 [2] CRAN (R 4.2.1)
# S4Vectors            * 0.36.2    2023-02-26 [2] Bioconductor
# ScaledMatrix           1.6.0     2022-11-01 [2] Bioconductor
# scales                 1.2.1     2022-08-20 [2] CRAN (R 4.2.1)
# scater               * 1.26.1    2022-11-13 [2] Bioconductor
# scuttle              * 1.8.4     2023-01-19 [2] Bioconductor
# sessioninfo            1.2.2     2021-12-06 [2] CRAN (R 4.2.1)
# SingleCellExperiment * 1.20.1    2023-03-17 [2] Bioconductor
# sparseMatrixStats      1.10.0    2022-11-01 [2] Bioconductor
# stringi                1.7.12    2023-01-11 [2] CRAN (R 4.2.2)
# stringr              * 1.5.0     2022-12-02 [2] CRAN (R 4.2.2)
# SummarizedExperiment * 1.28.0    2022-11-01 [2] Bioconductor
# tibble               * 3.2.1     2023-03-20 [2] CRAN (R 4.2.3)
# tidyr                * 1.3.0     2023-01-24 [2] CRAN (R 4.2.2)
# tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.2.1)
# tidyverse            * 2.0.0     2023-02-22 [2] CRAN (R 4.2.2)
# timechange             0.2.0     2023-01-11 [2] CRAN (R 4.2.2)
# tzdb                   0.3.0     2022-03-28 [2] CRAN (R 4.2.1)
# utf8                   1.2.3     2023-01-31 [2] CRAN (R 4.2.2)
# vctrs                  0.6.2     2023-04-19 [1] CRAN (R 4.2.3)
# vipor                  0.4.5     2017-03-22 [2] CRAN (R 4.2.1)
# viridis                0.6.2     2021-10-13 [2] CRAN (R 4.2.1)
# viridisLite            0.4.2     2023-05-02 [1] CRAN (R 4.2.3)
# withr                  2.5.0     2022-03-03 [2] CRAN (R 4.2.1)
# XVector                0.38.0    2022-11-01 [2] Bioconductor
# zlibbioc               1.44.0    2022-11-01 [2] Bioconductor
# 
# [1] /users/bsimbiat/R/4.2.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/library
