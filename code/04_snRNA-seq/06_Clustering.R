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
pdf(file = here("plots", "04_snRNA-seq", "06_Clustering", "TSNE_clusters_unfinalized.pdf"),
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
pdf(file = here("plots", "04_snRNA-seq", "06_Clustering", "UMAP_clusters_unfinalized.pdf"),
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
pdf(here("plots", "04_snRNA-seq", "06_Clustering",
         "Heatmap_Clusters_unfinalized_officialclusters.pdf"), height = 7, width = 7)
   
    lapply(colorbyGroup, function(n) {
      # adding Rand index (index of similary to ct_Erik)
      Rand <- signif(pairwiseRand(sce$ct_Erik, colData(sce)[, n], mode = "index"), 3)
      plot_cap <- paste("Rand Index Against Erik Cell Types:", Rand)
      
      linker <- linkClustersMatrix(sce$ct_Erik, colData(sce)[, n])
      
      Heatmap(linker,
              column_title = plot_cap,
              col = c("black", viridisLite::plasma(101)),
              name = "Corr"
      )
    }) 

dev.off()

############# COLORING BY DOUBLETS #############################################
pdf(file = here("plots", "04_snRNA-seq", "06_Clustering", "TSNE_clust_doubletscore_predrops.pdf"))


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

