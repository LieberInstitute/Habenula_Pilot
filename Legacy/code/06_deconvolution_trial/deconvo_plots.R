library(SummarizedExperiment)
library(RColorBrewer)
library(jaffelab)
library(here)
library(reshape2)
library(patchwork)
library(purrr)
library(tidyverse)
library(broom)
library(sessioninfo)
library(pheatmap)
library(SingleCellExperiment)

## Load colors and plotting functions
source(here("main_colors.R"))
source(here("deconvolution","big_little_boxplot.R"))
load(here("deconvolution","data","cell_colors.Rdata"), verbose = TRUE)

## Load MuSiC results & res data
load(here("deconvolution","data","est_prop_top5_long.Rdata"), verbose = TRUE)
load(here("deconvolution","data","est_prop_top5.Rdata"),verbose = TRUE)
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)

Dx_colors <- brewer.pal(3, "Dark2")
#### Boxplots ####
#walk2(est_prop_long, names(est_prop_long), function(x,y){
broad_boxplot(data = est_prop_long,
                            xvar = "cell_type", 
                            yvar = "prop",
                            fillvar =  "PrimaryDx",
                            colorvar = "ignore",
                            pallet = Dx_colors,
                            title = "MuSiC Proportions: Top5 markers",
                            subtitle = "Habenula_broad")
  ggsave(plot = blb, filename = paste0("plots/cellType_boxplots_",y,".png"), width = 15)
} )


dx_sex_colors <- mdd_Dx_colors_LD
names(dx_sex_colors) <- gsub("dark","F",gsub("light","M",names(mdd_Dx_colors_LD)))

walk2(est_prop_long, names(est_prop_long), function(x,y){
  
  x <- x %>% mutate(Dx_Sex = paste0(PrimaryDx, "_",Sex))
  
  blb <- big_little_boxplot(data = x,
                            xvar = "cell_type", 
                            yvar = "prop",
                            fillvar =  "Dx_Sex",
                            colorvar = "ignore",
                            pallet = dx_sex_colors,
                            title = "MuSiC Proptions: Top5 markers",
                            subtitle = y)
  ggsave(plot = blb, filename = paste0("plots/cellType_sex_boxplots_",y,".png"), width = 15)
} )


#### Experimental plots ####
bar <-est_prop_long[[1]]%>%
  left_join(est_prop_long[[1]] %>%
              filter(cell_type == 'Oligo') %>%
              select(sample, prop_1 = prop),
            by = 'sample') %>% 
  ggplot(aes(x = reorder(sample, desc(prop_1)), y = prop, fill = cell_type)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cell_colors)+
  coord_flip()

ggsave(bar, filename = "plots/bar_test.png", height = 10)

# ## try pheatmap
# my_anno_colors <- list(cellType = cell_colors,
#                        PrimaryDx = mdd_Dx_colors,
#                        Experiment = mdd_dataset_colors)
# 
# sample_anno <- as.data.frame(colData(rse_gene)) %>%
#   select(Experiment, Sex, PrimaryDx, AgeDeath)
# 
# walk2(est_prop, names(est_prop), function(est, n){
#   epw <- est$Est.prop.weighted
#   png(paste0("plots/MuSiC-heatmap_",n,".png"), height = 800, width = 580)
#   pheatmap(epw,
#            show_rownames = FALSE,
#            annotation_row = sample_anno,
#            annotation_colors = my_anno_colors,
#            main = paste(n,"cell type proptions"))
#   dev.off()
# })

#### Compare MuSiC output with sn data from same donors ####
load(here("deconvolution","data","sce.sacc_filtered.Rdata"), verbose = TRUE)
load(here("deconvolution","data","sce.amyg_filtered.Rdata"), verbose = TRUE)

donors <- unique(sce.sacc$donor)

input_cells <- list(sacc_broad = colData(sce.sacc)[,c("donor","cellType.Broad")],
                   amyg_broad = colData(sce.amyg)[,c("donor","cellType.Broad")],
                   sacc_specific = colData(sce.sacc)[,c("donor","cellType")],
                   amyg_specific = colData(sce.amyg)[,c("donor","cellType")])

input_prop <- map(input_cells, function(x){
  colnames(x)[2] <- "cell_type"
  
  x %>% as_tibble %>%
    group_by(donor, cell_type) %>%
    count() %>%
    group_by(donor) %>%
    mutate(input_prop = n/sum(n))
})

output_prop <- map2(est_prop_long, input_prop, function(est,input){
  est %>% select(donor = BrNum, cell_type, prop) %>%
    inner_join(input)
})

output_prop_df <- do.call("rbind", output_prop) %>%
  rownames_to_column(var = "dataset") %>%
  mutate(dataset = ss(dataset,"\\."))

# levels(output_prop_df$cell_type) <- levels(output_prop_df$cell_type)[order(levels(output_prop_df$cell_type))]

in_out_prop <- ggplot(output_prop_df, aes(input_prop, prop, color = cell_type, shape = donor))+
  geom_point() +
  scale_color_manual(values = cell_colors)+
  geom_abline()+
  facet_wrap(~dataset) +
  labs(y = "MuSiC estimated prop", x = "Input sn data prop")

ggsave(in_out_prop, filename = "plots/prop_input_vs_ouput.png")

# sgejobs::job_single('deconvo_plots', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript deconvo_plots.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# [1] "Reproducibility information:"
# [1] "2021-01-11 11:12:17 EST"
# user  system elapsed 
# 42.305   2.392  46.937 
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value                                      
# version  R version 4.0.3 Patched (2020-11-29 r79529)
# os       CentOS Linux 7 (Core)                      
# system   x86_64, linux-gnu                          
# ui       X11                                        
# language (EN)                                       
# collate  en_US.UTF-8                                
# ctype    en_US.UTF-8                                
# tz       US/Eastern                                 
# date     2021-01-11                                 
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version  date       lib source                                   
# assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.0.3)                           
# backports              1.2.0    2020-11-02 [1] CRAN (R 4.0.3)                           
# Biobase              * 2.50.0   2020-10-27 [2] Bioconductor                             
# BiocGenerics         * 0.36.0   2020-10-27 [2] Bioconductor                             
# bitops                 1.0-6    2013-08-17 [2] CRAN (R 4.0.3)                           
# broom                * 0.7.3    2020-12-16 [2] CRAN (R 4.0.3)                           
# cellranger             1.1.0    2016-07-27 [2] CRAN (R 4.0.3)                           
# cli                    2.2.0    2020-11-20 [1] CRAN (R 4.0.3)                           
# colorspace             2.0-0    2020-11-11 [2] CRAN (R 4.0.3)                           
# crayon                 1.3.4    2017-09-16 [2] CRAN (R 4.0.3)                           
# DBI                    1.1.0    2019-12-15 [2] CRAN (R 4.0.3)                           
# dbplyr                 2.0.0    2020-11-03 [2] CRAN (R 4.0.3)                           
# DelayedArray           0.16.0   2020-10-27 [2] Bioconductor                             
# digest                 0.6.27   2020-10-24 [1] CRAN (R 4.0.3)                           
# dplyr                * 1.0.2    2020-08-18 [1] CRAN (R 4.0.3)                           
# ellipsis               0.3.1    2020-05-15 [2] CRAN (R 4.0.3)                           
# fansi                  0.4.1    2020-01-08 [2] CRAN (R 4.0.3)                           
# farver                 2.0.3    2020-01-16 [2] CRAN (R 4.0.3)                           
# forcats              * 0.5.0    2020-03-01 [2] CRAN (R 4.0.3)                           
# fs                     1.5.0    2020-07-31 [1] CRAN (R 4.0.3)                           
# generics               0.1.0    2020-10-31 [2] CRAN (R 4.0.3)                           
# GenomeInfoDb         * 1.26.2   2020-12-08 [2] Bioconductor                             
# GenomeInfoDbData       1.2.4    2020-11-30 [2] Bioconductor                             
# GenomicRanges        * 1.42.0   2020-10-27 [2] Bioconductor                             
# ggplot2              * 3.3.3    2020-12-30 [2] CRAN (R 4.0.3)                           
# glue                   1.4.2    2020-08-27 [1] CRAN (R 4.0.3)                           
# googledrive            1.0.1    2020-05-05 [1] CRAN (R 4.0.3)                           
# gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.0.3)                           
# haven                  2.3.1    2020-06-01 [1] CRAN (R 4.0.3)                           
# here                 * 1.0.0    2020-11-15 [1] CRAN (R 4.0.3)                           
# hms                    0.5.3    2020-01-08 [2] CRAN (R 4.0.3)                           
# httr                   1.4.2    2020-07-20 [1] CRAN (R 4.0.3)                           
# IRanges              * 2.24.0   2020-10-27 [1] Bioconductor                             
# jaffelab             * 0.99.30  2020-11-02 [1] Github (LieberInstitute/jaffelab@42637ff)
# jsonlite               1.7.1    2020-09-07 [1] CRAN (R 4.0.3)                           
# labeling               0.4.2    2020-10-20 [2] CRAN (R 4.0.3)                           
# lattice                0.20-41  2020-04-02 [3] CRAN (R 4.0.3)                           
# lifecycle              0.2.0    2020-03-06 [1] CRAN (R 4.0.3)                           
# limma                  3.46.0   2020-10-27 [2] Bioconductor                             
# lubridate              1.7.9.2  2020-11-13 [1] CRAN (R 4.0.3)                           
# magrittr               2.0.1    2020-11-17 [2] CRAN (R 4.0.3)                           
# Matrix                 1.3-2    2021-01-06 [3] CRAN (R 4.0.3)                           
# MatrixGenerics       * 1.2.0    2020-10-27 [2] Bioconductor                             
# matrixStats          * 0.57.0   2020-09-25 [2] CRAN (R 4.0.3)                           
# modelr                 0.1.8    2020-05-19 [2] CRAN (R 4.0.3)                           
# munsell                0.5.0    2018-06-12 [2] CRAN (R 4.0.3)                           
# patchwork            * 1.1.0    2020-11-09 [1] CRAN (R 4.0.3)                           
# pheatmap             * 1.0.12   2019-01-04 [2] CRAN (R 4.0.3)                           
# pillar                 1.4.7    2020-11-20 [1] CRAN (R 4.0.3)                           
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.0.3)                           
# plyr                   1.8.6    2020-03-03 [2] CRAN (R 4.0.3)                           
# ps                     1.4.0    2020-10-07 [1] CRAN (R 4.0.3)                           
# purrr                * 0.3.4    2020-04-17 [1] CRAN (R 4.0.3)                           
# R6                     2.5.0    2020-10-28 [1] CRAN (R 4.0.3)                           
# rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 4.0.3)                           
# RColorBrewer         * 1.1-2    2014-12-07 [2] CRAN (R 4.0.3)                           
# Rcpp                   1.0.5    2020-07-06 [1] CRAN (R 4.0.3)                           
# RCurl                  1.98-1.2 2020-04-18 [2] CRAN (R 4.0.3)                           
# readr                * 1.4.0    2020-10-05 [2] CRAN (R 4.0.3)                           
# readxl                 1.3.1    2019-03-13 [2] CRAN (R 4.0.3)                           
# reprex                 0.3.0    2019-05-16 [2] CRAN (R 4.0.3)                           
# reshape2             * 1.4.4    2020-04-09 [2] CRAN (R 4.0.3)                           
# rlang                * 0.4.9    2020-11-26 [1] CRAN (R 4.0.3)                           
# rprojroot              2.0.2    2020-11-15 [2] CRAN (R 4.0.3)                           
# rstudioapi             0.13     2020-11-12 [2] CRAN (R 4.0.3)                           
# rvest                  0.3.6    2020-07-25 [2] CRAN (R 4.0.3)                           
# S4Vectors            * 0.28.1   2020-12-09 [2] Bioconductor                             
# scales                 1.1.1    2020-05-11 [2] CRAN (R 4.0.3)                           
# segmented              1.3-0    2020-10-27 [1] CRAN (R 4.0.3)                           
# sessioninfo          * 1.1.1    2018-11-05 [2] CRAN (R 4.0.3)                           
# SingleCellExperiment * 1.12.0   2020-10-27 [2] Bioconductor                             
# stringi                1.5.3    2020-09-09 [1] CRAN (R 4.0.3)                           
# stringr              * 1.4.0    2019-02-10 [2] CRAN (R 4.0.3)                           
# SummarizedExperiment * 1.20.0   2020-10-27 [2] Bioconductor                             
# tibble               * 3.0.4    2020-10-12 [1] CRAN (R 4.0.3)                           
# tidyr                * 1.1.2    2020-08-27 [2] CRAN (R 4.0.3)                           
# tidyselect             1.1.0    2020-05-11 [2] CRAN (R 4.0.3)                           
# tidyverse            * 1.3.0    2019-11-21 [1] CRAN (R 4.0.3)                           
# vctrs                  0.3.5    2020-11-17 [1] CRAN (R 4.0.3)                           
# withr                  2.3.0    2020-09-22 [1] CRAN (R 4.0.3)                           
# xml2                   1.3.2    2020-04-23 [2] CRAN (R 4.0.3)                           
# XVector                0.30.0   2020-10-27 [2] Bioconductor                             
# zlibbioc               1.36.0   2020-10-27 [2] Bioconductor                             
# 
# [1] /users/lhuuki/R/4.0.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/library
