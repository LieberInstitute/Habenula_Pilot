## May 15, 2023 - Bukola Ajanaku
# Probing deconvoluted data to find trends alongside certain phenotypes.
# qrsh -l mem_free=20G,h_vmem=20G

library(SummarizedExperiment)
library(here)
library(tidyverse)
library(ggplot2)

# loading deconvo data
load(file = here("processed-data", "99_paper_figs", "sce_objects", "prop_long.RDATA"),
     verbose = TRUE)