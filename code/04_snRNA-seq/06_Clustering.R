## 2/1/23 - Bukola Ajanaku
## Clustering post-harmony sce object (contains doublet scoring info) for clustering
## and future annotation.
## Based on:
# 1) https://github.com/LieberInstitute/Habenula_Bulk/blob/master/code/08_snRNA-seq_Erik/20210323_human_hb_neun.R
# 2) https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/code/03_build_sce/06_cluster.R
# 3) https://github.com/LieberInstitute/Habenula_Bulk/blob/master/code/09_snRNA-seq_re-processed/04_clustering.R
# 4) https://www.stephaniehicks.com/biocdemo/articles/Demo.html
# 5) OSCA workflows and multi-sample example sections: https://bioconductor.org/books/release/OSCA/book-contents.html
# qrsh -l mem_free=20G,h_vmem=20G

# Loading relevant libraries:
library("SingleCellExperiment")
library("jaffelab")
library("scater")
library("scran")
library("here")
library("sessioninfo")
library("ggplot2")
library("cowplot")
library("tidyverse")
library("RColorBrewer")
library("scales")
library("Polychrome")
library("testthat") # for Rand?
library("bluster") # for Rand
library("mbkmeans")
# library("pheatmap")
library("viridisLite")
library("ComplexHeatmap")

# Loading post harmony sce object 
load(here("processed-data", "04_snRNA-seq", "sce_objects", "sce_harmony_by_Samp.Rdata"))
sce <- sce_corrbySamp
rm(sce_corrbySamp)

##### Approach 1: Louvain Clustering (Steph Hicks) #############################
set.seed(777)
  firstWay <- buildSNNGraph(sce, k = 10, use.dimred = "HARMONY")
  lou <- igraph::cluster_louvain(firstWay)
  sce$louvain10 <- paste0("Louv_10_", lou$membership)
table(sce$louvain10)
  # Louv_10_1 Louv_10_10 Louv_10_11 Louv_10_12 Louv_10_13 Louv_10_14  Louv_10_2 
  # 1320        756        167       1754        542        350       1940 
  # Louv_10_3  Louv_10_4  Louv_10_5  Louv_10_6  Louv_10_7  Louv_10_8  Louv_10_9 
  # 1842       1607       2559       1315        745        514       1671 

set.seed(732)
  secondWay <- buildSNNGraph(sce, k = 20, use.dimred = "HARMONY")
  lou2 <- igraph::cluster_louvain(secondWay)
  sce$louvain20 <- paste0("Louv_20_", lou2$membership)
table(sce$louvain20)

  # Louv_20_1 Louv_20_10 Louv_20_11 Louv_20_12  Louv_20_2  Louv_20_3  Louv_20_4 
  # 1323        156       2167        544       2539       2022       1604 
  # Louv_20_5  Louv_20_6  Louv_20_7  Louv_20_8  Louv_20_9 
  # 1479       2039       1028       1431        750 

set.seed(702)
  thirdWay <- buildSNNGraph(sce, k = 50, use.dimred = "HARMONY")
  lou3 <- igraph::cluster_louvain(thirdWay)
  sce$louvain50 <- paste0("Louv_50_", lou3$membership)
table(sce$louvain50)
  # Louv_50_1 Louv_50_10  Louv_50_2  Louv_50_3  Louv_50_4  Louv_50_5  Louv_50_6 
  # 1381       2182       2384       2102       1668       2039       3068 
  # Louv_50_7  Louv_50_8  Louv_50_9 
  # 967        549        742 

##### Approach 2: WalkTrap Clustering (Erik) ###################################
set.seed(2565)
  g10 <- buildSNNGraph(sce, k=10, use.dimred = 'HARMONY')
  clust10 <- igraph::cluster_walktrap(g10)$membership
  colData(sce)$wT_10_Erik <-  paste0("10wTrap_", factor(clust10))
table(colData(sce)$wT_10_Erik)
    # 10wTrap_1 10wTrap_10 10wTrap_11 10wTrap_12 10wTrap_13 10wTrap_14 10wTrap_15 
    # 145       2231        188        152        201        325        164 
    # 10wTrap_16 10wTrap_17 10wTrap_18 10wTrap_19  10wTrap_2 10wTrap_20 10wTrap_21 
    # 145       5241        373         51       1796         86        217 
    # 10wTrap_22 10wTrap_23 10wTrap_24 10wTrap_25 10wTrap_26 10wTrap_27 10wTrap_28 
    # 276        280        266         62        134         65         65 
    # 10wTrap_29  10wTrap_3 10wTrap_30 10wTrap_31 10wTrap_32 10wTrap_33 10wTrap_34 
    # 85        477         51         83         38         39         25 
    # 10wTrap_35 10wTrap_36 10wTrap_37  10wTrap_4  10wTrap_5  10wTrap_6  10wTrap_7 
    # 17         18         17       1014        312         38       1833 
    # 10wTrap_8  10wTrap_9 
    # 177        395 

set.seed(1234)
  g20 <- buildSNNGraph(sce, k = 20, use.dimred = 'HARMONY')
  clust20 <- igraph::cluster_walktrap(g20)$membership
  sce$wT_20_Erik <- paste0("20wTrap_", factor(clust20))
table(sce$wT_20_Erik)
    # 20wTrap_1 20wTrap_10 20wTrap_11 20wTrap_12 20wTrap_13 20wTrap_14 20wTrap_15 
    # 2319        281        201       1062        613        189        110 
    # 20wTrap_16 20wTrap_17 20wTrap_18 20wTrap_19  20wTrap_2 20wTrap_20 20wTrap_21 
    # 283        184        227        289        723         34         75 
    # 20wTrap_22 20wTrap_23  20wTrap_3  20wTrap_4  20wTrap_5  20wTrap_6  20wTrap_7 
    # 75         30       1120        970        353       2695        335 
    # 20wTrap_8  20wTrap_9 
    # 4793        121 

set.seed(4321)
  g50 <- buildSNNGraph(sce, k=50, use.dimred = 'HARMONY')
  clust50 <- igraph::cluster_walktrap(g50)$membership
  colData(sce)$wT_50_Erik <- paste0("50wTrap_", factor(clust50))
table(colData(sce)$wT_50_Erik)
    # 50wTrap_1 50wTrap_10 50wTrap_11 50wTrap_12 50wTrap_13 50wTrap_14  50wTrap_2 
    # 2563       4470        679        301        125        257        510 
    # 50wTrap_3  50wTrap_4  50wTrap_5  50wTrap_6  50wTrap_7  50wTrap_8  50wTrap_9 
    # 1808        428       1427       3020        377        744        373 


##### Approach 3: Mini-batch k Means (Steph Hicks Example) #####################
  # set.seed(789)
  # k_list <- seq(5, 20)
  # km_res <- lapply(k_list, function(k) {
  #   mbkmeans(sce, clusters = k, 
  #            batch_size = 500,
  #            reduceMethod = "HARMONY",
  #            calc_wcss = TRUE)
  # })
  # wcss <- sapply(km_res, function(x) sum(x$WCSS_per_cluster))

# NOT PROCEEDING WITH MBK

############# Pre-PLOTTING #####################################################
## NOTE: 
## 1) Getting rid of walkTrap 5 because we do not need 100+ grups 

# Fixing NAs in ct_Erik
sce$ct_Erik <- as.character(sce$ct_Erik)
sce$ct_Erik[is.na(sce$ct_Erik)] <- "Bukola"
    # Astro       Bukola         Endo        LHb.1        LHb.2        LHb.3 
    # 495          887           78         1344          744          614 
    # LHb.4        LHb.5        LHb.6        MHb.1        MHb.2        Micro 
    # 494          302           72          487          197          317 
    # Neuron.Ambig      Oligo.1      Oligo.2      OPC.1        OPC.2  Thal.GABA.1 
    # 5                  1680         426          220          745         2632 
    # Thal.GABA.2  Thal.GABA.3      Thal.MD      Thal.PF     Thal.PVT 
    # 4269               80          232          220          542 

## Fixing factors for sce$ct_Erik
sce$ct_Erik <- factor(sce$ct_Erik, levels = c("Astro", "Endo", "Micro", "Oligo.1", "Oligo.2", "OPC.1", 
                        "OPC.2", "LHb.1", "LHb.2", "LHb.3", "LHb.4", "LHb.5", "LHb.6", 
                        "MHb.1", "MHb.2", "Neuron.Ambig", "Thal.GABA.1", 
                        "Thal.GABA.2", "Thal.GABA.3", "Thal.MD", "Thal.PF", "Thal.PVT", 
                        "Bukola"))

# Grouping Erik Cell Types into Particular Groups
sce$groupErik <- NA

colData(sce)[sce$ct_Erik %in% c("Astro", "Endo", "Micro", "Oligo.1", "Oligo.2", "OPC.1", "OPC.2"), 
             "groupErik"] <- "Glia_and_Endo"
colData(sce)[sce$ct_Erik %in% c("LHb.1", "LHb.2", "LHb.3", "LHb.4", "LHb.5", "LHb.6", 
                                "MHb.1", "MHb.2"), "groupErik"] <- "Hb_Neurons"
colData(sce)[sce$ct_Erik %in% c("Thal.GABA.1", "Thal.GABA.2", "Thal.GABA.3", "Thal.MD",
                                "Thal.PF", "Thal.PVT"), "groupErik"] <- "Thal_Neurons"
colData(sce)[sce$ct_Erik == "Neuron.Ambig", "groupErik"] <- "Erik_Neuron.Ambig"
colData(sce)[sce$ct_Erik == "Bukola", "groupErik"] <- "Bukola_New"

any(is.na(sce$groupErik))
# FALSE

# Creating list to color by 
colorbyGroup <- c(colnames(colData(sce)[grepl("wT_", colnames(colData(sce)))]), 
                  colnames(colData(sce)[grepl("louvain", colnames(colData(sce)))]))
  # [1] "wT_10_Erik" "wT_20_Erik" "wT_50_Erik" "louvain10"  "louvain20" 
  # [6] "louvain50"

# creating a color palette
dark2 = c("#A6BDD7", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", 
          "#A6761D", "#666666")

############# PLOTTING #########################################################
############# REGULAR TSNE & UMAPS #############################################
# Plotting harmonized TSNE with colorbyGroup
pdf(file = here("plots", "04_snRNA-seq", "06_Clustering", "TSNE_clusters_unfinalized_officialclusters.pdf"),
   height = 6,  width = 13)

lapply(colorbyGroup, function(n) {
  # figuring our color scale
  if (length(table(colData(sce)[,n])) <= 8){
    
    colorer <- dark2
  
    } else if (length(table(colData(sce)[,n])) > 8){
    colourCount = length(table(colData(sce)[,n]))
    getPalette = colorRampPalette(brewer.pal(9, "Set1"))
    
    colorer <-  getPalette(colourCount)
    }
  
  # adding Rand index (index of similary on scale over 1)
  Rand <- signif(pairwiseRand(sce$ct_Erik, colData(sce)[, n], mode = "index"), 3)
  plot_cap <- paste("Rand Index:", Rand)
    
  # plotting full plot
  regplot <- ggcells(sce, mapping = aes(x = TSNE.1, y = TSNE.2, colour = .data[[n]])) +
    geom_point(size = 0.9, alpha = 0.5)  +
    coord_equal() +
    labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2") + 
   ggtitle(paste("Total Number of Groups =", length(table(colData(sce)[,n])))) +
    scale_colour_manual(values = colorer) + theme(legend.position = "None") + 
    theme_classic()
  
  # creating split by original cell type Erik for right hand side
  faceted <- ggcells(sce, mapping = aes(x = TSNE.1, y = TSNE.2, colour = .data[[n]])) +
      geom_point(size = 0.5, alpha = 0.5)  +
      coord_equal() +
      labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2") + facet_wrap(~ ct_Erik) +
      guides(colour = guide_legend(override.aes = list(size = 3))) + 
      scale_colour_manual(values = colorer) +  
      theme_classic() 

  ggdraw(add_sub(plot_grid(regplot, faceted, nrow = 1), plot_cap))
  
})

dev.off()


# Plotting harmonized UMAPS with colorbyGroup 
pdf(file = here("plots", "04_snRNA-seq", "06_Clustering", "UMAP_clusters_unfinalized_officialclusters.pdf"),
    height = 7, width = 12)

lapply(colorbyGroup, function(n) {
  if (length(table(colData(sce)[,n])) <= 8){
    
    colorer <- dark2
    
  } else if (length(table(colData(sce)[,n])) > 8){
    colourCount = length(table(colData(sce)[,n]))
    getPalette = colorRampPalette(brewer.pal(9, "Set1"))
    
    colorer <-  getPalette(colourCount)
  }
  
  # adding Rand index (index of similary on scale over 1)
  Rand <- signif(pairwiseRand(sce$ct_Erik, colData(sce)[, n], mode = "index"), 3)
  plot_cap <- paste("Rand Index:", Rand)
  
  # plotting full plot
  regplot <- ggcells(sce, mapping = aes(x = UMAP.1, y = UMAP.2, colour = .data[[n]])) +
    geom_point(size = 0.9, alpha = 0.5)  +
    coord_equal() +
    labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2") + 
    ggtitle(paste("Total Number of Groups =", length(table(colData(sce)[,n])))) +
    scale_colour_manual(values = colorer) + theme(legend.position = "None") +
    theme_classic()
  
  # creating split by original cell type Erik for right hand side
  faceted <- ggcells(sce, mapping = aes(x = UMAP.1, y = UMAP.2, colour = .data[[n]])) +
    geom_point(size = 0.4, alpha = 0.5)  +
    coord_equal() +
    labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2") + facet_wrap(~ ct_Erik) +
    guides(colour = guide_legend(override.aes = list(size = 3))) + 
    scale_colour_manual(values = colorer) +  theme_classic()
  
  plot_grid(regplot, faceted, nrow = 1)
  
  ggdraw(add_sub(plot_grid(regplot, faceted, nrow = 1), plot_cap))
  
})

dev.off()

############# HEATMAPS #########################################################
# With Dr. Torres' coloring scheme
pdf(here("plots", "04_snRNA-seq", "06_Clustering",
         "Heatmap_Clusters_unfinalized_officialclusters.pdf"), height = 7, width = 7)
   
    lapply(colorbyGroup, function(n) {
      # adding Rand index (index of similary to ct_Erik)
      Rand <- signif(pairwiseRand(sce$ct_Erik, colData(sce)[, n], mode = "index"), 3)
      plot_cap <- paste("Rand Index Against Erik Cell Types:", Rand)
      numClusts <- paste("Number of Clusters =", length(table(colData(sce)[,n])))
      
      linker <- linkClustersMatrix(sce$ct_Erik, colData(sce)[, n])
      
      Heatmap(linker,
              column_title = plot_cap,
              col = c("black", viridisLite::plasma(101)),
              name = "Corr",
              row_title = numClusts
      )
    }) 

dev.off()

############# COLORING BY DOUBLETS #############################################
pdf(file = here("plots", "04_snRNA-seq", "06_Clustering", "TSNE_clust_doubletscore_predrops_officialclusters.pdf"))

  # plotting full plot
  regplot <- ggcells(sce, mapping = aes(x = TSNE.1, y = TSNE.2, colour = .data[["doubletScore"]])) +
    geom_point(size = 0.9, alpha = 0.5)  +
    coord_equal() +
    labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2") + 
    ggtitle("TSNE by Doublet Score") +
    theme(legend.position = "Right") + 
    theme_classic() 
  
  # creating split by original cell type Erik for right hand side
  faceted <-ggcells(sce, mapping = aes(x = TSNE.1, y = TSNE.2, colour = .data[["doubletScore"]])) +
    geom_point(size = 0.9, alpha = 0.5)  +
    coord_equal() +
    labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2") + 
    ggtitle("TSNE by Doublet Score") +
    theme(legend.position = "Right") + 
    theme_classic() +  
    facet_wrap(~ ct_Erik) 
   
  regplot
  faceted 
 
dev.off()

### SAVING OBJECT(S) ###########################################################
save(sce, file = here("processed-data", "04_snRNA-seq", "sce_objects", 
                      "sce_post_clustering.Rdata"))


## Reproducibility information
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
# beachmat               2.14.2    2023-04-07 [2] Bioconductor
# beeswarm               0.4.0     2021-06-01 [2] CRAN (R 4.2.1)
# benchmarkme            1.0.8     2022-06-12 [2] CRAN (R 4.2.1)
# benchmarkmeData        1.0.4     2020-04-23 [2] CRAN (R 4.2.1)
# Biobase              * 2.58.0    2022-11-01 [2] Bioconductor
# BiocGenerics         * 0.44.0    2022-11-01 [2] Bioconductor
# BiocNeighbors          1.16.0    2022-11-01 [2] Bioconductor
# BiocParallel           1.32.6    2023-03-17 [2] Bioconductor
# BiocSingular           1.14.0    2022-11-01 [2] Bioconductor
# bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.2.1)
# bluster              * 1.8.0     2022-11-01 [2] Bioconductor
# brio                   1.1.3     2021-11-30 [2] CRAN (R 4.2.1)
# circlize               0.4.15    2022-05-10 [2] CRAN (R 4.2.1)
# cli                    3.6.1     2023-03-23 [2] CRAN (R 4.2.3)
# clue                   0.3-64    2023-01-31 [2] CRAN (R 4.2.2)
# cluster                2.1.4     2022-08-22 [3] CRAN (R 4.2.3)
# ClusterR               1.3.0     2023-01-21 [2] CRAN (R 4.2.2)
# codetools              0.2-19    2023-02-01 [3] CRAN (R 4.2.3)
# colorout             * 1.2-2     2023-02-09 [1] Github (jalvesaq/colorout@79931fd)
# colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.2.2)
# ComplexHeatmap       * 2.14.0    2022-11-01 [2] Bioconductor
# cowplot              * 1.1.1     2020-12-30 [2] CRAN (R 4.2.1)
# crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.2.1)
# DelayedArray           0.24.0    2022-11-01 [2] Bioconductor
# DelayedMatrixStats     1.20.0    2022-11-01 [2] Bioconductor
# digest                 0.6.31    2022-12-11 [2] CRAN (R 4.2.2)
# doParallel             1.0.17    2022-02-07 [2] CRAN (R 4.2.1)
# dplyr                * 1.1.1     2023-03-22 [2] CRAN (R 4.2.3)
# dqrng                  0.3.0     2021-05-01 [2] CRAN (R 4.2.1)
# edgeR                  3.40.2    2023-01-19 [2] Bioconductor
# fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.2.2)
# forcats              * 1.0.0     2023-01-29 [2] CRAN (R 4.2.2)
# foreach                1.5.2     2022-02-02 [2] CRAN (R 4.2.1)
# fs                     1.6.1     2023-02-06 [2] CRAN (R 4.2.2)
# gargle                 1.3.0     2023-01-30 [2] CRAN (R 4.2.2)
# generics               0.1.3     2022-07-05 [2] CRAN (R 4.2.1)
# GenomeInfoDb         * 1.34.9    2023-02-02 [2] Bioconductor
# GenomeInfoDbData       1.2.9     2022-09-29 [2] Bioconductor
# GenomicRanges        * 1.50.2    2022-12-16 [2] Bioconductor
# GetoptLong             1.0.5     2020-12-15 [2] CRAN (R 4.2.1)
# ggbeeswarm             0.7.1     2022-12-16 [2] CRAN (R 4.2.2)
# ggplot2              * 3.4.2     2023-04-03 [2] CRAN (R 4.2.3)
# ggrepel                0.9.3     2023-02-03 [2] CRAN (R 4.2.2)
# GlobalOptions          0.1.2     2020-06-10 [2] CRAN (R 4.2.1)
# glue                   1.6.2     2022-02-24 [2] CRAN (R 4.2.1)
# gmp                    0.7-1     2023-02-07 [2] CRAN (R 4.2.2)
# googledrive            2.1.0     2023-03-22 [2] CRAN (R 4.2.3)
# gridExtra              2.3       2017-09-09 [2] CRAN (R 4.2.1)
# gtable                 0.3.3     2023-03-21 [2] CRAN (R 4.2.3)
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.2.1)
# hms                    1.1.3     2023-03-21 [2] CRAN (R 4.2.3)
# httr                   1.4.5     2023-02-24 [2] CRAN (R 4.2.2)
# igraph                 1.4.2     2023-04-07 [2] CRAN (R 4.2.3)
# IRanges              * 2.32.0    2022-11-01 [2] Bioconductor
# irlba                  2.3.5.1   2022-10-03 [2] CRAN (R 4.2.1)
# iterators              1.0.14    2022-02-05 [2] CRAN (R 4.2.1)
# jaffelab             * 0.99.32   2023-02-09 [1] Github (LieberInstitute/jaffelab@7b7afe3)
# lattice                0.20-45   2021-09-22 [3] CRAN (R 4.2.3)
# lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.2.1)
# limma                  3.54.2    2023-02-28 [2] Bioconductor
# locfit                 1.5-9.7   2023-01-02 [2] CRAN (R 4.2.2)
# lubridate            * 1.9.2     2023-02-10 [2] CRAN (R 4.2.2)
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.2.1)
# MASS                   7.3-58.2  2023-01-23 [3] CRAN (R 4.2.3)
# Matrix                 1.5-4     2023-04-04 [2] CRAN (R 4.2.3)
# MatrixGenerics       * 1.10.0    2022-11-01 [2] Bioconductor
# matrixStats          * 1.0.0     2023-06-02 [1] CRAN (R 4.2.3)
# mbkmeans             * 1.14.0    2022-11-01 [2] Bioconductor
# metapod                1.6.0     2022-11-01 [2] Bioconductor
# munsell                0.5.0     2018-06-12 [2] CRAN (R 4.2.1)
# nlme                   3.1-162   2023-01-31 [2] CRAN (R 4.2.2)
# pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.2.3)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.2.1)
# png                    0.1-8     2022-11-29 [2] CRAN (R 4.2.2)
# Polychrome           * 1.5.1     2022-05-03 [1] CRAN (R 4.2.2)
# purrr                * 1.0.1     2023-01-10 [2] CRAN (R 4.2.2)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.2.1)
# rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 4.2.2)
# RColorBrewer         * 1.1-3     2022-04-03 [2] CRAN (R 4.2.1)
# Rcpp                   1.0.10    2023-01-22 [2] CRAN (R 4.2.2)
# RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.2.3)
# readr                * 2.1.4     2023-02-10 [2] CRAN (R 4.2.2)
# rjson                  0.2.21    2022-01-09 [2] CRAN (R 4.2.1)
# rlang                  1.1.1     2023-04-28 [1] CRAN (R 4.2.3)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.2.1)
# rsvd                   1.0.5     2021-04-16 [2] CRAN (R 4.2.1)
# S4Vectors            * 0.36.2    2023-02-26 [2] Bioconductor
# ScaledMatrix           1.6.0     2022-11-01 [2] Bioconductor
# scales               * 1.2.1     2022-08-20 [2] CRAN (R 4.2.1)
# scater               * 1.26.1    2022-11-13 [2] Bioconductor
# scatterplot3d          0.3-44    2023-05-05 [1] CRAN (R 4.2.3)
# scran                * 1.26.2    2023-01-19 [2] Bioconductor
# scuttle              * 1.8.4     2023-01-19 [2] Bioconductor
# segmented              1.6-4     2023-04-13 [1] CRAN (R 4.2.3)
# sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.2.1)
# shape                  1.4.6     2021-05-19 [2] CRAN (R 4.2.1)
# SingleCellExperiment * 1.20.1    2023-03-17 [2] Bioconductor
# sparseMatrixStats      1.10.0    2022-11-01 [2] Bioconductor
# statmod                1.5.0     2023-01-06 [2] CRAN (R 4.2.2)
# stringi                1.7.12    2023-01-11 [2] CRAN (R 4.2.2)
# stringr              * 1.5.0     2022-12-02 [2] CRAN (R 4.2.2)
# SummarizedExperiment * 1.28.0    2022-11-01 [2] Bioconductor
# testthat             * 3.1.7     2023-03-12 [2] CRAN (R 4.2.3)
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
# viridisLite          * 0.4.2     2023-05-02 [1] CRAN (R 4.2.3)
# withr                  2.5.0     2022-03-03 [2] CRAN (R 4.2.1)
# XVector                0.38.0    2022-11-01 [2] Bioconductor
# zlibbioc               1.44.0    2022-11-01 [2] Bioconductor
# 
# [1] /users/bsimbiat/R/4.2.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/library

