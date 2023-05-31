## May 2, 2023 - Bukola Ajanaku
# Using PCs, UMAPs, and TSNEs for supplement fig for thal vs habenula.
# qrsh -l mem_free=20G,h_vmem=20G

library("SingleCellExperiment")
library("here")
library("tidyverse")
library("scater")
library("ggplot2")
library("cowplot")

# Loading sce object
load(here("processed-data", "99_paper_figs", "sce_objects", 
          "official_final_sce.RDATA"), verbose = TRUE)
# sce

table(sce$final_Annotations)
# Astrocyte       Endo Excit.Thal Inhib.Thal  LHb.1      LHb.2      LHb.3 
# 538         38       1800       7612        201        266        134 
# LHb.4      LHb.5      LHb.6      LHb.7      MHb.1      MHb.2      MHb.3 
# 477         83         39       1014        152        540         18 
# Microglia      Oligo        OPC 
# 145       2178       1796 

# sourcing official color palette 
source(file = here("code", "99_paper_figs", "source_colors.R"))
# bulk_colors and sn_colors

# creating plot directory
plot_dir <- here("plots", "99_paper_figs", "07_Hb_vs_Thal_Plots")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

##### PLOTTING  ################################################################
# post-harmony PC2 vs PC3
pdf(here(plot_dir, "PC2_vs_PC3_splitbyClusters.pdf"),
    width = 15, height = 8)
plot1 <- plotReducedDim(sce, dimred = "PCA", ncomponents = 2:3) +
  geom_point(aes(color = sce$final_Annotations)) +
  scale_colour_manual(values = sn_colors) +
  theme(legend.position = "none")

plot2 <- plotReducedDim(sce, dimred = "PCA", ncomponents = 2:3) +
  geom_point(aes(color = sce$final_Annotations)) +
  scale_colour_manual(values = sn_colors) +
  facet_wrap(~ sce$final_Annotations) +
  guides(color = guide_legend(title="Cell Type"))

plot_grid(plot1, plot2)
dev.off()


##### FOR ONE DRIVE PIECES #####################################################
# pdf 
pdf(here(plot_dir, "forOneDrive","sfigu_PC2vsPC3_HabvsThalModality.pdf"),
    width = 8, height = 4.5)
    plot1 <- plotReducedDim(sce, dimred = "PCA", ncomponents = 2:3) +
      geom_point(aes(color = sce$final_Annotations)) +
      scale_colour_manual(values = sn_colors) +
      theme(legend.position = "none")
    
    plot2 <- plotReducedDim(sce, dimred = "PCA", ncomponents = 2:3) +
      geom_point(aes(color = sce$final_Annotations)) +
      scale_colour_manual(values = sn_colors) +
      facet_wrap(~ sce$final_Annotations) +
      guides(color = guide_legend(title="Cell Type")) + 
      theme(axis.text = element_text(size = 6)) 

  plot_grid(plot1, plot2,
            rel_widths = c(.4, .6))
dev.off()

# png
png(here(plot_dir, "forOneDrive","sfigu_PC2vsPC3_HabvsThalModality.png"),
    width = 8, height = 4.5, units = "in", res = 1200)
plot1 <- plotReducedDim(sce, dimred = "PCA", ncomponents = 2:3) +
  geom_point(aes(color = sce$final_Annotations)) +
  scale_colour_manual(values = sn_colors) +
  theme(legend.position = "none")

plot2 <- plotReducedDim(sce, dimred = "PCA", ncomponents = 2:3) +
  geom_point(aes(color = sce$final_Annotations)) +
  scale_colour_manual(values = sn_colors) +
  facet_wrap(~ sce$final_Annotations) +
  guides(color = guide_legend(title="Cell Type")) + 
  theme(axis.text = element_text(size = 6)) 

plot_grid(plot1, plot2,
          rel_widths = c(.4, .6))
dev.off()









# Done.