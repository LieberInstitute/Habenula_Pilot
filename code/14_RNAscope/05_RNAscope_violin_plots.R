
library("tidyverse")
library("GGally")
library("here")
library("sessioninfo")
library("patchwork")

## set up
plot_dir <- here("plots", "14_RNAscope", "02_MHb_HALO")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
