# 00_Color_Scheme_ct.R 
# 2/28/23 - Bukola Ajanaku (advised and mentored by Louise Huuki-Myers)
## This creates a custom color scheme based on colors grabbed from 
## https://medialab.github.io/iwanthue/ . This operates as a function that returns 
## a list of HEX colors that can change lengths based on number of colors needed.
## Max amount for now will be 20.

library("here")
library("jaffelab")
library("RColorBrewer")
library("sessioninfo")

color_pallete <- c("#5b679a",
  "#c2cf1a",
  "#d779ff",
  "#009421",
  "#78319a",
  "#00af80",
  "#ff4b7f",
  "#0091e1",
  "#d45a00",
  "#919fff",
  "#ffb050",
  "#faa0ff",
  "#525b00",
  "#ff64b8",
  "#806200",
  "#ffa5c7",
  "#91332d",
  "#ff886b",
  "#8b3160",
  "#b8766a")

grabColors <- function(num, start = 1){
  # num is the total number of colors needed
  # start defaults at the first color but can be re-written to start elsewhere
  return(color_pallete[start:(start + num - 1)])
}


## How Louise wanted me to do this:

# cellTypecolors_9 <-
#   c(ExcNeuron = "#274ea6",
#     InhibThalMed = "#6ed94f",
#     OPC = "#ff67e3",
#     Astro = "#296f00",
#     Microglia = "#cb0029",
#     ThalMedPVT = "#009770",
#     Oligo = "#ff815f",
#     ExcThalMed = "#d794b4",
#     LHbN = "#ffb848")

preview_colors <- function(cell_colors) {
  par(las = 2) # make label text perpendicular to axis
  par(mar = c(5, 8, 4, 2)) # increase y-axis margin.
  barplot(rep(1, length(cell_colors)),
          col = cell_colors,
          horiz = TRUE,
          axes = FALSE,
          names.arg = names(cell_colors)
  )
}

 png(here("plots", "cell_colors", "cellTypecolors_9.png"), height = 800)
   preview_colors(cellTypecolors_9)
# preview_colors(grabColors(9))
 dev.off()

# # saving colors
# save(cellTypecolors_9, file = here("processed-data", "cell_type_colors.Rdata"))

# Reproducibility Information:
print("Reproducibility information:")
options(width = 120)
session_info()
  
# ─ Session info ─────────────────────────────────────────────────────────────────────────────────────────────────────────────
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
# ─ Packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# Biobase                2.58.0    2022-11-01 [2] Bioconductor
# BiocGenerics           0.44.0    2022-11-01 [2] Bioconductor
# bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.2.1)
# cli                    3.6.1     2023-03-23 [2] CRAN (R 4.2.3)
# colorout             * 1.2-2     2023-02-09 [1] Github (jalvesaq/colorout@79931fd)
# DelayedArray           0.24.0    2022-11-01 [2] Bioconductor
# dplyr                  1.1.1     2023-03-22 [2] CRAN (R 4.2.3)
# fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.2.2)
# fs                     1.6.1     2023-02-06 [2] CRAN (R 4.2.2)
# gargle                 1.3.0     2023-01-30 [2] CRAN (R 4.2.2)
# generics               0.1.3     2022-07-05 [2] CRAN (R 4.2.1)
# GenomeInfoDb           1.34.9    2023-02-02 [2] Bioconductor
# GenomeInfoDbData       1.2.9     2022-09-29 [2] Bioconductor
# GenomicRanges          1.50.2    2022-12-16 [2] Bioconductor
# glue                   1.6.2     2022-02-24 [2] CRAN (R 4.2.1)
# googledrive            2.1.0     2023-03-22 [2] CRAN (R 4.2.3)
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.2.1)
# IRanges                2.32.0    2022-11-01 [2] Bioconductor
# jaffelab             * 0.99.32   2023-02-09 [1] Github (LieberInstitute/jaffelab@7b7afe3)
# lattice                0.20-45   2021-09-22 [3] CRAN (R 4.2.3)
# lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.2.1)
# limma                  3.54.2    2023-02-28 [2] Bioconductor
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.2.1)
# MASS                   7.3-58.2  2023-01-23 [3] CRAN (R 4.2.3)
# Matrix                 1.5-4     2023-04-04 [2] CRAN (R 4.2.3)
# MatrixGenerics         1.10.0    2022-11-01 [2] Bioconductor
# matrixStats            1.0.0     2023-06-02 [1] CRAN (R 4.2.3)
# nlme                   3.1-162   2023-01-31 [2] CRAN (R 4.2.2)
# pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.2.3)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.2.1)
# purrr                  1.0.1     2023-01-10 [2] CRAN (R 4.2.2)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.2.1)
# rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 4.2.2)
# RColorBrewer         * 1.1-3     2022-04-03 [2] CRAN (R 4.2.1)
# RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.2.3)
# rlang                  1.1.1     2023-04-28 [1] CRAN (R 4.2.3)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.2.1)
# S4Vectors              0.36.2    2023-02-26 [2] Bioconductor
# segmented              1.6-4     2023-04-13 [1] CRAN (R 4.2.3)
# sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.2.1)
# SummarizedExperiment   1.28.0    2022-11-01 [2] Bioconductor
# tibble                 3.2.1     2023-03-20 [2] CRAN (R 4.2.3)
# tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.2.1)
# utf8                   1.2.3     2023-01-31 [2] CRAN (R 4.2.2)
# vctrs                  0.6.2     2023-04-19 [1] CRAN (R 4.2.3)
# XVector                0.38.0    2022-11-01 [2] Bioconductor
# zlibbioc               1.44.0    2022-11-01 [2] Bioconductor
# 
# [1] /users/bsimbiat/R/4.2.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/library
# 
