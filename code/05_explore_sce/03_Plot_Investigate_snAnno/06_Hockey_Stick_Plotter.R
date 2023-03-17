## March 10, 2023 - Bukola Ajanaku
# Making hockey stick plots for mean expression data on snAnno2 and snAnno3 
# qrsh -l mem_free=20G,h_vmem=20G

library("SingleCellExperiment")
library("here")
library("tidyverse")
library("sessioninfo")
library("DeconvoBuddies")
library("ggplot2")

# loading regular mean ratio and 1vAll data 
load(here("processed-data", "05_explore_sce", "mean_ratio_for_snAnno3_from_05_Updated_Annotations_meanExpression.Rdata"),
          verbose = TRUE)

# setting up plot_dir
plot_dir <- here("plots", "05_explore_sce", "03_Plot_Investigate_snAnno",
                 "06_Hockey_Stick_Plotter")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

# sourcing code from DLPFC Project (by Louise Huuki) 
source(here("code", "04_snRNA-seq", "sourcing", "custom_plotExpression.R"))

# actually runs and plots markers (by Louise Huuki) 
source(here("code", "04_snRNA-seq", "sourcing", "my_plotMarkers.R"))

# sourcing my grabColors(), max colors 20 (by Bukola Ajanaku) 
source(here("code", "04_snRNA-seq", "sourcing", "color_Scheme_CT.R"))

# adding color group
snAnno_marker_stats_new_3_combined$Top25 <- "No"
snAnno_marker_stats_new_3_combined[which(snAnno_marker_stats_new_3_combined$rank_ratio <= 25),
                                          "Top25"] <- "Yes"

# plotting
# this is the hockey stick plot for the split up snAnno
pdf(here("plots", "05_explore_sce", "03_Plot_Investigate_snAnno",
         "06_Hockey_Stick_Plotter", "ratio_vs_stdFC-snAnno3_freescales_2.pdf"))
   ggplot(snAnno_marker_stats_new_3_combined, aes(ratio, std.logFC)) +
     geom_point(size = 0.5, aes(colour = Top25)) +  
  facet_wrap(~cellType.target, scales = "free") 
dev.off()


# labs(x = "mean(target logcount)/mean(highest non-target logcount)") +

# ahhh, I get it now. We cannot do hockey stick plots on snAnno. Time to move on
# to next folder because this is now getting into the bulk deconvo terrain.
