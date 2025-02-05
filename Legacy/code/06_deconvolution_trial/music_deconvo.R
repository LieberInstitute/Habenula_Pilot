
library(SummarizedExperiment)
library(SingleCellExperiment)
library(jaffelab)
library(here)
library(MuSiC)
library(Biobase)
library(xbioc)
library(purrr)
library(tidyverse)
library(reshape2)
library(sessioninfo)

#### Load Data ####
## Load rse data
load(here("count_data","rse_gene_Roche_Habenula_qcAndAnnotated_n69.Rdata"), 
     verbose = TRUE)
rownames(rse_gene) <- rowData(rse_gene)$ensemblID
rse_gene<-rse_gene[unique(rownames(rse_gene)),]

## sce Data
load(here("Deconvolution","data","sce.all.n12_filtered.Rdata"), verbose = TRUE)
sce.all.n12<-sce.all.n12[unique(rownames(sce.all.n12)),]
sce.all.n12<-sce.all.n12[,!sce.all.n12$cellType.Broad=='Ambig']
colData(sce.all.n12)$cellType.Broad<-droplevels(colData(sce.all.n12)$cellType.Broad)

## marker data
load(here("Deconvolution","data","marker_stats.Rdata"), verbose = TRUE)
top_n <- 5

marker_genes <- marker_stats %>%
                      arrange(ratio) %>%
                      filter(rank_ratio <= top_n) %>%
                      arrange(cellType.target, rank_ratio) %>% 
                      pull(gene)

map_int(marker_genes, length)
map_int(marker_genes, ~sum(.x %in% rownames(rse_gene)))


exp <- list(hb = list(bulk = rse_gene[,rse_gene$Brain.Region == "Habenula"],
                        sce = sce.all.n12) 
  
)
map(exp, ~map_int(.x, ncol))
# $sacc
# bulk  sce 
# 551 7004 
# 
# $amyg
# bulk  sce 
# 540 6582 
exp_set <- map(exp, ~list(bulk = ExpressionSet(assayData = assays(.x$bulk)$counts),
                          sce = ExpressionSet(assayData = as.matrix(assays(.x$sce)$counts),
                                              phenoData=AnnotatedDataFrame(
                                                as.data.frame(colData(.x$sce)
                                              )
                          )
               )))
      

#### estimate cell type props ####
#ct <- list(broad = "cellType.Broad", specific = "cellType")

#est_prop <- map2(marker_genes, names(marker_genes), function(mg, n){
est_prop <- music_prop(bulk.eset = exp_set$hb$bulk,
             sc.eset = exp_set$hb$sce,
             clusters = phenoData(exp_set$hb$sce)[[24]],
             samples = 'donor',
             markers = marker_genes)
#return(est_prop)
#} 
#)
colMeans(est_prop$Est.prop.weighted)
#map(est_prop, ~round(colMeans(.x$Est.prop.weighted),3))
# $sacc_broad
# Oligo Micro Astro Inhib   OPC Excit 
# 0.156 0.052 0.552 0.211 0.000 0.028 
# 
# $amyg_broad
# Inhib Oligo Astro Excit Micro   OPC 
# 0.192 0.169 0.474 0.085 0.072 0.008 
# 
# $sacc_specific
# Oligo   Micro   Astro Inhib.2 Inhib.1     OPC Excit.3 Excit.1 Excit.2 Excit.4 
# 0.220   0.080   0.595   0.011   0.035   0.003   0.017   0.038   0.001   0.000 
# 
# $amyg_specific
# Inhib.2   Oligo   Astro Excit.1   Micro     OPC Inhib.3 Inhib.5 Inhib.1 
# 0.023   0.296   0.348   0.051   0.189   0.003   0.001   0.053   0.035

save(est_prop, file= here("Deconvolution", "data",
                          paste0("est_prop_top", top_n,".Rdata")))

pd <- as.data.frame(colData(rse_gene))

pd2 <- pd  %>%
  select(BrNum, RNum, Sex, PrimaryDx, ERCCsumLogErr, grep("snpPC",colnames(pd)))%>%
  rownames_to_column("sample")
pd2$sample<-as.integer(pd2$sample)  

est_prop_long <- melt(est_prop$Est.prop.weighted) %>%
                       rename(sample = Var1, cell_type = Var2, prop = value) %>% 
                       arrange(cell_type) %>%
                       left_join(pd2, by = "sample")
est_prop_long$cell_type<-factor(est_prop_long$cell_type)
#levels(est_prop_long$cell_type)<-levels(est_prop_long$cell_type)[order(levels(est_prop_long$cell_type))]

## order celltypes
#est_prop_long <- map(est_prop_long, function(x){
  #alphabetize
#  l_alpha <- levels(x$cell_type)[order(levels(x$cell_type))]
#  i_inhib <- grep("Inhib",l_alpha)
#  i_excit <- grep("Excit",l_alpha)
#  l_order <- c(l_alpha[!1:length(l_alpha) %in% c(i_inhib, i_excit)],l_alpha[c(i_excit, i_inhib)])
#  levels(x$cell_type) <- l_order
#  return(x)
#})

save(est_prop_long, file = here("Deconvolution","data",paste0("est_prop_top", top_n,"_long.Rdata")))

# sgejobs::job_single('music_deconvo', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript music_deconvo.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# [1] "Reproducibility information:"
# [1] "2020-12-21 12:15:50 EST"
# user  system elapsed 
# 69.009   6.005  77.247 
# �?? Session info �??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??
# setting  value                                      
# version  R version 4.0.3 Patched (2020-11-29 r79529)
# os       CentOS Linux 7 (Core)                      
# system   x86_64, linux-gnu                          
# ui       X11                                        
# language (EN)                                       
# collate  en_US.UTF-8                                
# ctype    en_US.UTF-8                                
# tz       US/Eastern                                 
# date     2020-12-21                                 
# 
# �?? Packages �??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??�??
# package              * version    date       lib source                                   
# AnnotationDbi        * 1.52.0     2020-10-27 [2] Bioconductor                             
# assertthat             0.2.1      2019-03-21 [2] CRAN (R 4.0.3)                           
# backports              1.2.0      2020-11-02 [1] CRAN (R 4.0.3)                           
# Biobase              * 2.50.0     2020-10-27 [2] Bioconductor                             
# BiocGenerics         * 0.36.0     2020-10-27 [2] Bioconductor                             
# BiocManager            1.30.10    2019-11-16 [2] CRAN (R 4.0.3)                           
# bit                    4.0.4      2020-08-04 [2] CRAN (R 4.0.3)                           
# bit64                  4.0.5      2020-08-30 [2] CRAN (R 4.0.3)                           
# bitops                 1.0-6      2013-08-17 [2] CRAN (R 4.0.3)                           
# blob                   1.2.1      2020-01-20 [2] CRAN (R 4.0.3)                           
# broom                  0.7.3      2020-12-16 [2] CRAN (R 4.0.3)                           
# cellranger             1.1.0      2016-07-27 [2] CRAN (R 4.0.3)                           
# checkmate              2.0.0      2020-02-06 [2] CRAN (R 4.0.3)                           
# cli                    2.2.0      2020-11-20 [1] CRAN (R 4.0.3)                           
# coda                   0.19-4     2020-09-30 [2] CRAN (R 4.0.3)                           
# codetools              0.2-18     2020-11-04 [3] CRAN (R 4.0.3)                           
# colorspace             2.0-0      2020-11-11 [2] CRAN (R 4.0.3)                           
# conquer                1.0.2      2020-08-27 [2] CRAN (R 4.0.3)                           
# crayon                 1.3.4      2017-09-16 [2] CRAN (R 4.0.3)                           
# DBI                    1.1.0      2019-12-15 [2] CRAN (R 4.0.3)                           
# dbplyr                 2.0.0      2020-11-03 [2] CRAN (R 4.0.3)                           
# DelayedArray           0.16.0     2020-10-27 [2] Bioconductor                             
# digest                 0.6.27     2020-10-24 [1] CRAN (R 4.0.3)                           
# dplyr                * 1.0.2      2020-08-18 [1] CRAN (R 4.0.3)                           
# ellipsis               0.3.1      2020-05-15 [2] CRAN (R 4.0.3)                           
# fansi                  0.4.1      2020-01-08 [2] CRAN (R 4.0.3)                           
# forcats              * 0.5.0      2020-03-01 [2] CRAN (R 4.0.3)                           
# fs                     1.5.0      2020-07-31 [1] CRAN (R 4.0.3)                           
# generics               0.1.0      2020-10-31 [2] CRAN (R 4.0.3)                           
# GenomeInfoDb         * 1.26.2     2020-12-08 [2] Bioconductor                             
# GenomeInfoDbData       1.2.4      2020-11-30 [2] Bioconductor                             
# GenomicRanges        * 1.42.0     2020-10-27 [2] Bioconductor                             
# ggplot2              * 3.3.2      2020-06-19 [2] CRAN (R 4.0.3)                           
# glue                   1.4.2      2020-08-27 [1] CRAN (R 4.0.3)                           
# googledrive            1.0.1      2020-05-05 [1] CRAN (R 4.0.3)                           
# gtable                 0.3.0      2019-03-25 [2] CRAN (R 4.0.3)                           
# haven                  2.3.1      2020-06-01 [1] CRAN (R 4.0.3)                           
# here                 * 1.0.0      2020-11-15 [1] CRAN (R 4.0.3)                           
# hms                    0.5.3      2020-01-08 [2] CRAN (R 4.0.3)                           
# httr                   1.4.2      2020-07-20 [1] CRAN (R 4.0.3)                           
# IRanges              * 2.24.0     2020-10-27 [1] Bioconductor                             
# jaffelab             * 0.99.30    2020-11-02 [1] Github (LieberInstitute/jaffelab@42637ff)
# jsonlite               1.7.1      2020-09-07 [1] CRAN (R 4.0.3)                           
# lattice                0.20-41    2020-04-02 [3] CRAN (R 4.0.3)                           
# lifecycle              0.2.0      2020-03-06 [1] CRAN (R 4.0.3)                           
# limma                  3.46.0     2020-10-27 [2] Bioconductor                             
# lubridate              1.7.9.2    2020-11-13 [1] CRAN (R 4.0.3)                           
# magrittr               2.0.1      2020-11-17 [2] CRAN (R 4.0.3)                           
# MASS                   7.3-53     2020-09-09 [3] CRAN (R 4.0.3)                           
# Matrix                 1.2-18     2019-11-27 [3] CRAN (R 4.0.3)                           
# MatrixGenerics       * 1.2.0      2020-10-27 [2] Bioconductor                             
# MatrixModels           0.4-1      2015-08-22 [2] CRAN (R 4.0.3)                           
# matrixStats          * 0.57.0     2020-09-25 [2] CRAN (R 4.0.3)                           
# mcmc                   0.9-7      2020-03-21 [2] CRAN (R 4.0.3)                           
# MCMCpack               1.4-9      2020-08-02 [2] CRAN (R 4.0.3)                           
# memoise                1.1.0      2017-04-21 [2] CRAN (R 4.0.3)                           
# modelr                 0.1.8      2020-05-19 [2] CRAN (R 4.0.3)                           
# munsell                0.5.0      2018-06-12 [2] CRAN (R 4.0.3)                           
# MuSiC                * 0.1.1      2020-11-02 [1] Github (xuranw/MuSiC@01e51ba)            
# nnls                 * 1.4        2012-03-19 [1] CRAN (R 4.0.3)                           
# pillar                 1.4.7      2020-11-20 [1] CRAN (R 4.0.3)                           
# pkgconfig              2.0.3      2019-09-22 [2] CRAN (R 4.0.3)                           
# pkgmaker               0.32.2.900 2020-11-02 [1] Github (renozao/pkgmaker@51b7207)        
# plyr                   1.8.6      2020-03-03 [2] CRAN (R 4.0.3)                           
# ps                     1.4.0      2020-10-07 [1] CRAN (R 4.0.3)                           
# purrr                * 0.3.4      2020-04-17 [1] CRAN (R 4.0.3)                           
# quantreg               5.75       2020-10-27 [2] CRAN (R 4.0.3)                           
# R6                     2.5.0      2020-10-28 [1] CRAN (R 4.0.3)                           
# rafalib              * 1.0.0      2015-08-09 [1] CRAN (R 4.0.3)                           
# RColorBrewer           1.1-2      2014-12-07 [2] CRAN (R 4.0.3)                           
# Rcpp                   1.0.5      2020-07-06 [1] CRAN (R 4.0.3)                           
# RCurl                  1.98-1.2   2020-04-18 [2] CRAN (R 4.0.3)                           
# readr                * 1.4.0      2020-10-05 [2] CRAN (R 4.0.3)                           
# readxl                 1.3.1      2019-03-13 [2] CRAN (R 4.0.3)                           
# registry               0.5-1      2019-03-05 [2] CRAN (R 4.0.3)                           
# reprex                 0.3.0      2019-05-16 [2] CRAN (R 4.0.3)                           
# reshape2             * 1.4.4      2020-04-09 [2] CRAN (R 4.0.3)                           
# rlang                  0.4.9      2020-11-26 [1] CRAN (R 4.0.3)                           
# rprojroot              2.0.2      2020-11-15 [2] CRAN (R 4.0.3)                           
# RSQLite                2.2.1      2020-09-30 [2] CRAN (R 4.0.3)                           
# rstudioapi             0.13       2020-11-12 [2] CRAN (R 4.0.3)                           
# rvest                  0.3.6      2020-07-25 [2] CRAN (R 4.0.3)                           
# S4Vectors            * 0.28.1     2020-12-09 [2] Bioconductor                             
# scales                 1.1.1      2020-05-11 [2] CRAN (R 4.0.3)                           
# segmented              1.3-0      2020-10-27 [1] CRAN (R 4.0.3)                           
# sessioninfo          * 1.1.1      2018-11-05 [2] CRAN (R 4.0.3)                           
# SingleCellExperiment * 1.12.0     2020-10-27 [2] Bioconductor                             
# SparseM                1.78       2019-12-13 [2] CRAN (R 4.0.3)                           
# stringi                1.5.3      2020-09-09 [1] CRAN (R 4.0.3)                           
# stringr              * 1.4.0      2019-02-10 [2] CRAN (R 4.0.3)                           
# SummarizedExperiment * 1.20.0     2020-10-27 [2] Bioconductor                             
# tibble               * 3.0.4      2020-10-12 [1] CRAN (R 4.0.3)                           
# tidyr                * 1.1.2      2020-08-27 [2] CRAN (R 4.0.3)                           
# tidyselect             1.1.0      2020-05-11 [2] CRAN (R 4.0.3)                           
# tidyverse            * 1.3.0      2019-11-21 [1] CRAN (R 4.0.3)                           
# vctrs                  0.3.5      2020-11-17 [1] CRAN (R 4.0.3)                           
# withr                  2.3.0      2020-09-22 [1] CRAN (R 4.0.3)                           
# xbioc                * 0.1.19     2020-11-02 [1] Github (renozao/xbioc@1354168)           
# xml2                   1.3.2      2020-04-23 [2] CRAN (R 4.0.3)                           
# xtable                 1.8-4      2019-04-21 [2] CRAN (R 4.0.3)                           
# XVector                0.30.0     2020-10-27 [2] Bioconductor                             
# zlibbioc               1.36.0     2020-10-27 [2] Bioconductor                             
# 
# [1] /users/lhuuki/R/4.0.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/library
# **** Job ends ****
#   Mon Dec 21 12:15:51 EST 2020

