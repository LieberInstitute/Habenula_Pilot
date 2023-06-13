## May 15, 2023 - Bukola Ajanaku
# Probing deconvoluted data to find trends alongside certain phenotypes.
# qrsh -l mem_free=20G,h_vmem=20G

library(SummarizedExperiment)
library(here)
library(tidyverse)
library(ggplot2)
library(jaffelab)
library(ggrepel)
library(cowplot)
library(ggthemes)
library(viridis)
library(RColorBrewer)

# loading deconvo data
load(file = here("processed-data", "99_paper_figs", "sce_objects", "prop_long.RDATA"),
     verbose = TRUE)

head(prop_long, 10)
    # Sample cellType     prop BrNum PrimaryDx factor_CT  Hb_sum Br_Order prop_perc
    # <chr>  <chr>       <dbl> <chr> <chr>     <fct>       <dbl> <fct>        <dbl>
    #   1 R18424 Astrocyte  0      Br55… Control   Astrocyte 0       Br5572        0   
    # 2 R18424 Endo       0      Br55… Control   Endo      0       Br5572        0   
    # 3 R18424 Microglia  0      Br55… Control   Microglia 0       Br5572        0   
    # 4 R18424 Oligo      0      Br55… Control   Oligo     0       Br5572        0   
    # 5 R18424 OPC        0      Br55… Control   OPC       0       Br5572        0   
    # 6 R18424 Inhib.Thal 1      Br55… Control   Inhib.Th… 0       Br5572      100   
    # 7 R18424 Excit.Thal 0      Br55… Control   Excit.Th… 0       Br5572        0   
    # 8 R18424 MHb        0      Br55… Control   MHb       0       Br5572        0   
    # 9 R18424 LHb        0      Br55… Control   LHb       0       Br5572        0   
    # 10 R18358 Astrocyte  0.0335 Br14… Schizo    Astrocyte 0.00897 Br1427        3.35

# loading final rse object 
load(here("processed-data", "02_bulk_qc", "count_data_bukola", 
          "rse_gene_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata"), verbose = TRUE)
# rse_gene

# creating plotting directory
plot_dir <- here("plots", "99_paper_figs", "11_bulk_Investigation")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

# adding Thalamus proportions and then dividing Habenula over Thalamus
sum_Prop <- prop_long |>
  filter(cellType %in% c("Excit.Thal", "Inhib.Thal")) |>
  group_by(BrNum) |>
  summarize(Thal_sum = sum(prop)) 

glia_Prop <- prop_long |>
  filter(cellType %in% c("Astrocyte", "Endo", "Microglia", "Oligo", "OPC")) |>
  group_by(BrNum) |>
  summarize(Glia_sum = sum(prop)) 

prop_long <- left_join(prop_long, sum_Prop) |>
  arrange(Br_Order) |>
  mutate(Hb_over_Thal = (Hb_sum / Thal_sum) ) |>
  mutate(prop_perc = NULL) 

prop_long <- left_join(prop_long, glia_Prop) |>
  arrange(Br_Order)
    # %>% print(width = Inf)

# making data frame for ease of merge with rse ColData
comp_invest <- as.data.frame(prop_long[,c("Sample", "Hb_sum", "Thal_sum", "Hb_over_Thal",
                                          "Glia_sum", "BrNum")])
comp_invest$RNum = comp_invest$Sample


new_rse_colData <- as_tibble(colData(rse_gene))
new_rse_colData <- new_rse_colData |>
                    left_join(comp_invest[match(unique(comp_invest$RNum), comp_invest$RNum),],
                              copy = TRUE
                              )

new_rse_colData <- DataFrame(new_rse_colData)

# setting new colData 
colData(rse_gene) <- new_rse_colData

############## GRABBING PC VALUES ##############################################
# grabbing PC values 
pca <- prcomp(t(assays(rse_gene)$logcounts))
  
  ## % of the variance explained by each PC
pca_vars <- getPcaVars(pca)
pca_vars_labs<- paste0("PC", seq(along = pca_vars), ": ",
                       pca_vars, "% Var Expl")

  ## Joining PC and sample info
pca_data<-cbind(pca$x, colData(rse_gene))
  
  ## Add samples' phenotypes
pca_data<-as.data.frame(pca_data)
  
##### PLOTTING #################################################################
phenoMets <- c("AgeDeath", "PrimaryDx", "Flowcell", "mitoRate",
               "RIN")

cont_vars <- c("AgeDeath", "mitoRate", "RIN") 
disc_vars <- c("PrimaryDx", "Flowcell")

# adding no glia problem sanples 
pca_data$Notes <- "Normal"

# only thalamus 
pca_data[pca_data$BrNum %in% c("Br5572"), ]$Notes <- "All.Thal"

# no glia 
pca_data[pca_data$BrNum %in% c("Br5558", "Br2015", "Br8050", "Br1682", 
                               "Br6197"), ]$Notes <- "NoGlia"


# creating function for plotting
plot_pcs <- function(xer, yer, color_by){ 
  pos = position_jitter(width = 0.3, height = 0.3, seed = 1)
  
  plot_list = list()
  c = 1
  
  for(i in color_by) { 
    
    if(i %in% disc_vars){ 

      
  
        plot_list[[c]] <- ggplot(pca_data, 
                                 aes_string(x = xer,
                                            y = yer)) +
                                  geom_point() +
                                  geom_jitter(
                                    aes_string(x = xer,
                                               y = yer,
                                               color = as.factor(pca_data[, i])),
                                    position = pos,
                                    size = 2) +
                                  geom_text_repel(position = pos,
                                                  color = "darkgrey",
                                                  max.overlaps = 7,
                                                  size = 4,
                                                  aes_string(label = "BrNum"),
                                                  data = pca_data[pca_data$Notes == "Normal",]
                                  )  + 
                                  geom_text_repel(aes_string(label = "BrNum"),
                                                  color = "red", 
                                                  data = pca_data[pca_data$Notes == "NoGlia",],
                                                  position = pos,
                                                  size = 4) +
                                  geom_text_repel(aes_string(label = "BrNum"),
                                                  color = "black", 
                                                  data = pca_data[pca_data$Notes == "All.Thal",],
                                                  position = pos,
                                                  size = 4) +
                                  ggtitle(paste(xer, "vs", yer, "Colored by", i)) + 
                                  guides(color =guide_legend(title= i)) +
                                  scale_color_brewer(palette = "Dark2")
                          
    } else if(i %in% cont_vars){ 
      
      plot_list[[c]] <- ggplot(pca_data, 
                               aes_string(x = xer,
                                          y = yer)) +
                              geom_point() + 
                              geom_jitter(
                                position = pos,
                                size = 2,
                                aes_string(color = pca_data[, i])) +
                              geom_text_repel(position = pos,
                                              color = "darkgrey",
                                              max.overlaps = 7,
                                              size = 4,
                                              aes_string(label = "BrNum"),
                                              data = pca_data[pca_data$Notes == "Normal",]
                              )  + 
                              geom_text_repel(aes_string(label = "BrNum"),
                                              color = "red", 
                                              data = pca_data[pca_data$Notes == "NoGlia",],
                                              position = pos,
                                              size = 4) +
                              geom_text_repel(aes_string(label = "BrNum"),
                                              color = "black", 
                                              data = pca_data[pca_data$Notes == "All.Thal",],
                                              position = pos,
                                              size = 4) +
                              ggtitle(paste(xer, "vs", yer, "Colored by", i)) + 
                              guides(color = guide_legend(title= i)) +
                              scale_color_distiller(palette = "YlGnBu")
                      
      }
    
    
      c = c + 1
  }
  
  return(plot_list) 
} 

## Plotting PC vs PC
pdf(file = here(plot_dir, "PCvsPC", "PC1_VS_PC2_Investigation.pdf")) 
  plot_pcs(x = "PC1", y = "PC2", color_by = phenoMets)
dev.off()

pdf(file = here(plot_dir, "PCvsPC", "PC2_VS_PC3_Investigation.pdf")) 
  plot_pcs(x = "PC2", y = "PC3", color_by = phenoMets)
dev.off()

# creating similar function that plots metrics in a specific manner 
plot_mets <- function(xer, qc_met, color_by, metric_title){ 
  
  pos = position_jitter(seed = 777)
  
  plot_list = list()
  c = 1
  for(i in color_by) { 
    
    if(i %in% disc_vars){ 

      plot_list[[c]] <- ggplot(pca_data, 
                               aes_string(x = xer,
                                          y = qc_met)) +
        geom_point() +
        geom_jitter(
          aes_string(x = xer,
                     y = qc_met,
                     color = as.factor(pca_data[, i])),
                     position = pos,
                     size = 2) +
        geom_text_repel(position = pos,
                        color = "darkgrey",
                        max.overlaps = 7,
                        size = 4,
                        aes_string(label = "BrNum"),
                        data = pca_data[pca_data$Notes == "Normal",]
        )  + 
        geom_text_repel(aes_string(label = "BrNum"),
                        color = "red", 
                        data = pca_data[pca_data$Notes == "NoGlia",],
                        position = pos,
                        size = 4) +
        geom_text_repel(aes_string(label = "BrNum"),
                        color = "black", 
                        data = pca_data[pca_data$Notes == "All.Thal",],
                        position = pos,
                        size = 4) +
        ggtitle(paste(metric_title, "by", i)) + 
        guides(color =guide_legend(title= i)) +
        scale_color_brewer(palette = "Dark2") +
        ylab(metric_title)
      
    } else if(i %in% cont_vars){ 
      
      plot_list[[c]] <- ggplot(pca_data, 
                               aes_string(x = xer,
                                          y = qc_met)) +
        geom_point() + 
        geom_jitter(
          position = pos,
          size = 2,
          aes_string(color = pca_data[, i])) +
        geom_text_repel(position = pos,
                        color = "darkgrey",
                        max.overlaps = 7,
                        size = 4,
                        aes_string(label = "BrNum"),
                        data = pca_data[pca_data$Notes == "Normal",]
        )  + 
        geom_text_repel(aes_string(label = "BrNum"),
                        color = "red", 
                        data = pca_data[pca_data$Notes == "NoGlia",],
                        position = pos,
                        size = 4) +
        geom_text_repel(aes_string(label = "BrNum"),
                        color = "black", 
                        data = pca_data[pca_data$Notes == "All.Thal",],
                        position = pos,
                        size = 4) +
        ggtitle(paste(metric_title, "by", i)) + 
        guides(color = guide_legend(title= i)) +
        scale_color_distiller(palette = "YlGnBu") + 
        ylab(metric_title)
      
    }
    
    
    c = c + 1
  }
  
  return(plot_list) 
} 


# metrics against PC values 
# PC1
pdf(file = here(plot_dir, "Metric_Against_PC", "Hb_over_Thal_VS_PC1_Investigation.pdf")) 
  plot_mets(xer = "PC1", qc_met = "Hb_over_Thal", color_by = phenoMets, 
            metric_title = "Hb/Thal")
dev.off()

pdf(file = here(plot_dir, "Metric_Against_PC", "Hb_sum_VS_PC1_Investigation.pdf")) 
  plot_mets(xer = "PC1", qc_met = "Hb_sum", color_by = phenoMets, 
            metric_title = "Hb Sum Ratio")
dev.off()
pdf(file = here(plot_dir, "Metric_Against_PC", "Thal_sum_VS_PC1_Investigation.pdf")) 
  plot_mets(xer = "PC1", qc_met = "Thal_sum", color_by = phenoMets, 
            metric_title = "Thal Sum Ratio")
dev.off()

# PC2
pdf(file = here(plot_dir, "Metric_Against_PC", "Hb_over_Thal_VS_PC2_Investigation.pdf")) 
  plot_mets(xer = "PC2", qc_met = "Hb_over_Thal", color_by = phenoMets, 
            metric_title = "Hb/Thal")
dev.off()

pdf(file = here(plot_dir, "Metric_Against_PC", "Hb_sum_VS_PC2_Investigation.pdf")) 
  plot_mets(xer = "PC2", qc_met = "Hb_sum", color_by = phenoMets, 
            metric_title = "Hb Sum Ratio")
dev.off()

pdf(file = here(plot_dir, "Metric_Against_PC", "Thal_sum_VS_PC2_Investigation.pdf")) 
  plot_mets(xer = "PC2", qc_met = "Thal_sum", color_by = phenoMets, 
            metric_title = "Thal Sum Ratio")
dev.off()

# PC3:
pdf(file = here(plot_dir, "Metric_Against_PC", "Hb_over_Thal_VS_PC3_Investigation.pdf")) 
  plot_mets(xer = "PC3", qc_met = "Hb_over_Thal", color_by = phenoMets, 
            metric_title = "Hb/Thal")
dev.off()

pdf(file = here(plot_dir, "Metric_Against_PC", "Hb_sum_VS_PC3_Investigation.pdf")) 
  plot_mets(xer = "PC3", qc_met = "Hb_sum", color_by = phenoMets, 
            metric_title = "Hb Sum Ratio")
dev.off()

pdf(file = here(plot_dir, "Metric_Against_PC", "Thal_sum_VS_PC3_Investigation.pdf")) 
  plot_mets(xer = "PC3", qc_met = "Thal_sum", color_by = phenoMets, 
            metric_title = "Thal Sum Ratio")
dev.off()

## Metrics Against Glia Sum 
pdf(file = here(plot_dir, "Metric_Against_Glia", "Hb_sum_VS_Glia_sum.pdf")) 
  plot_mets(xer = "Glia_sum", qc_met = "Hb_sum", color_by = phenoMets, 
            metric_title = "Hb Sum Ratio")
dev.off()

pdf(file = here(plot_dir, "Metric_Against_Glia", "Thal_sum_VS_Glia_sum.pdf")) 
plot_mets(xer = "Glia_sum", qc_met = "Thal_sum", color_by = phenoMets, 
          metric_title = "Thal Sum Ratio")
dev.off()

pdf(file = here(plot_dir, "Metric_Against_Glia", "Hb_over_Thal_VS_Glia_sum.pdf")) 
plot_mets(xer = "Glia_sum", qc_met = "Hb_over_Thal", color_by = phenoMets, 
          metric_title = "Hb/Thal")
dev.off()

# pc vs glia
pdf(file = here(plot_dir, "Metric_Against_Glia", "PC3_VS_Glia_sum.pdf")) 
plot_mets(xer = "PC3", qc_met = "Glia_sum", color_by = phenoMets, 
          metric_title = "Glia Sum Prop")
dev.off()

pdf(file = here(plot_dir, "Metric_Against_Glia", "PC2_VS_Glia_sum.pdf")) 
plot_mets(xer = "PC2", qc_met = "Glia_sum", color_by = phenoMets, 
          metric_title = "Glia Sum Prop")
dev.off()

pdf(file = here(plot_dir, "Metric_Against_Glia", "PC1_VS_Glia_sum.pdf")) 
plot_mets(xer = "PC1", qc_met = "Glia_sum", color_by = phenoMets, 
          metric_title = "Glia Sum Prop")
dev.off()

pdf(file = here(plot_dir, "Metric_Against_Glia", "PC4_VS_Glia_sum.pdf")) 
plot_mets(xer = "PC4", qc_met = "Glia_sum", color_by = phenoMets, 
          metric_title = "Glia Sum Prop")
dev.off()

pdf(file = here(plot_dir, "Metric_Against_Glia", "PC5_VS_Glia_sum.pdf")) 
plot_mets(xer = "PC5", qc_met = "Glia_sum", color_by = phenoMets, 
          metric_title = "Glia Sum Prop")
dev.off()
# 

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
# Biobase              * 2.58.0    2022-11-01 [2] Bioconductor
# BiocGenerics         * 0.44.0    2022-11-01 [2] Bioconductor
# bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.2.1)
# cli                    3.6.1     2023-03-23 [2] CRAN (R 4.2.3)
# colorout             * 1.2-2     2023-02-09 [1] Github (jalvesaq/colorout@79931fd)
# colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.2.2)
# cowplot              * 1.1.1     2020-12-30 [2] CRAN (R 4.2.1)
# DelayedArray           0.24.0    2022-11-01 [2] Bioconductor
# dplyr                * 1.1.1     2023-03-22 [2] CRAN (R 4.2.3)
# fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.2.2)
# forcats              * 1.0.0     2023-01-29 [2] CRAN (R 4.2.2)
# fs                     1.6.1     2023-02-06 [2] CRAN (R 4.2.2)
# gargle                 1.3.0     2023-01-30 [2] CRAN (R 4.2.2)
# generics               0.1.3     2022-07-05 [2] CRAN (R 4.2.1)
# GenomeInfoDb         * 1.34.9    2023-02-02 [2] Bioconductor
# GenomeInfoDbData       1.2.9     2022-09-29 [2] Bioconductor
# GenomicRanges        * 1.50.2    2022-12-16 [2] Bioconductor
# ggplot2              * 3.4.2     2023-04-03 [2] CRAN (R 4.2.3)
# ggrepel              * 0.9.3     2023-02-03 [2] CRAN (R 4.2.2)
# ggthemes             * 4.2.4     2021-01-20 [1] CRAN (R 4.2.3)
# glue                   1.6.2     2022-02-24 [2] CRAN (R 4.2.1)
# googledrive            2.1.0     2023-03-22 [2] CRAN (R 4.2.3)
# gridExtra              2.3       2017-09-09 [2] CRAN (R 4.2.1)
# gtable                 0.3.3     2023-03-21 [2] CRAN (R 4.2.3)
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.2.1)
# hms                    1.1.3     2023-03-21 [2] CRAN (R 4.2.3)
# IRanges              * 2.32.0    2022-11-01 [2] Bioconductor
# jaffelab             * 0.99.32   2023-02-09 [1] Github (LieberInstitute/jaffelab@7b7afe3)
# lattice                0.20-45   2021-09-22 [3] CRAN (R 4.2.3)
# lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.2.1)
# limma                  3.54.2    2023-02-28 [2] Bioconductor
# lubridate            * 1.9.2     2023-02-10 [2] CRAN (R 4.2.2)
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.2.1)
# MASS                   7.3-58.2  2023-01-23 [3] CRAN (R 4.2.3)
# Matrix                 1.5-4     2023-04-04 [2] CRAN (R 4.2.3)
# MatrixGenerics       * 1.10.0    2022-11-01 [2] Bioconductor
# matrixStats          * 1.0.0     2023-06-02 [1] CRAN (R 4.2.3)
# munsell                0.5.0     2018-06-12 [2] CRAN (R 4.2.1)
# nlme                   3.1-162   2023-01-31 [2] CRAN (R 4.2.2)
# pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.2.3)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.2.1)
# purrr                * 1.0.1     2023-01-10 [2] CRAN (R 4.2.2)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.2.1)
# rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 4.2.2)
# RColorBrewer         * 1.1-3     2022-04-03 [2] CRAN (R 4.2.1)
# Rcpp                   1.0.10    2023-01-22 [2] CRAN (R 4.2.2)
# RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.2.3)
# readr                * 2.1.4     2023-02-10 [2] CRAN (R 4.2.2)
# rlang                  1.1.1     2023-04-28 [1] CRAN (R 4.2.3)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.2.1)
# S4Vectors            * 0.36.2    2023-02-26 [2] Bioconductor
# scales                 1.2.1     2022-08-20 [2] CRAN (R 4.2.1)
# segmented              1.6-4     2023-04-13 [1] CRAN (R 4.2.3)
# sessioninfo            1.2.2     2021-12-06 [2] CRAN (R 4.2.1)
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
# viridis              * 0.6.2     2021-10-13 [2] CRAN (R 4.2.1)
# viridisLite          * 0.4.2     2023-05-02 [1] CRAN (R 4.2.3)
# withr                  2.5.0     2022-03-03 [2] CRAN (R 4.2.1)
# XVector                0.38.0    2022-11-01 [2] Bioconductor
# zlibbioc               1.44.0    2022-11-01 [2] Bioconductor
# 
# [1] /users/bsimbiat/R/4.2.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/library
