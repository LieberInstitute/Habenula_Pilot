## March 10, 2023 - Bukola Ajanaku
# Making hockey stick plots for mean expression data on snAnno2 and snAnno3 
# qrsh -l mem_free=20G,h_vmem=20G

library("SingleCellExperiment")
library("here")
library("tidyverse")
library("sessioninfo")
library("DeconvoBuddies")

# loading regular sce object
load(here("processed-data", "05_explore_sce", 
          "mean_ratio_for_snAnno3_from_05_Updated_Annotations_meanExpression.Rdata",
          verbose = TRUE))


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

# 

