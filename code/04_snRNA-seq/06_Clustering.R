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


# Meeting with Dr. Torres (2/13):
# 1) change colors for heatmaps
# 2) plot gene markers for group 50 as well
# 3) get right dimensions for graphs,  replot all
# 4) prep for meeting tomorrrow at 3

# Loading post harmony sce object 
load(here("processed-data", "04_snRNA-seq", "sce_objects", "sce_harmony_by_Samp.Rdata"))
sce <- sce_corrbySamp
rm(sce_corrbySamp)

##### Approach 1: Louvain Clustering (Steph Hicks) #############################

set.seed(777)

firstWay <- buildSNNGraph(sce, k = 10, use.dimred = "HARMONY")
lou <- igraph::cluster_louvain(firstWay)
sce$louvain <- paste0("Louvain", lou$membership)
table(sce$louvain)
# Louvain1 Louvain10 Louvain11 Louvain12 Louvain13 Louvain14  Louvain2  Louvain3 
# 1328       756       178      1798       543       347      2254      1832 
# Louvain4  Louvain5  Louvain6  Louvain7  Louvain8  Louvain9 
# 1116      1855       620      1643      1089      1723


##### Approach 2: WalkTrap Clustering (Erik) ###################################
set.seed(1234)
g20 <- buildSNNGraph(sce, k = 20, use.dimred = 'HARMONY')
clust20 <- igraph::cluster_walktrap(g20)$membership
sce$k_20_Erik <- paste0("20wTrap", factor(clust20))
table(sce$k_20_Erik)
    # 20wTrap1 20wTrap10 20wTrap11 20wTrap12 20wTrap13 20wTrap14 20wTrap15 20wTrap16 
    # 1299       127       214      4747       296       192       151       347 
    # 20wTrap17 20wTrap18 20wTrap19  20wTrap2 20wTrap20 20wTrap21 20wTrap22  20wTrap3 
    # 211        34        75      1909        86       175        29      1867 
    # 20wTrap4  20wTrap5  20wTrap6  20wTrap7  20wTrap8  20wTrap9 
    # 569       450       359      2742       344       859 

set.seed(4321)
g50 <- buildSNNGraph(sce, k=50, use.dimred = 'HARMONY')
clust50 <- igraph::cluster_walktrap(g50)$membership
colData(sce)$k_50_Erik <- paste0("50wTrap", factor(clust50))
table(colData(sce)$k_50_Erik)
    # 50wTrap1 50wTrap10 50wTrap11 50wTrap12 50wTrap13 50wTrap14  50wTrap2  50wTrap3 
    # 110      4320       725       350       275       274       538       958 
    # 50wTrap4  50wTrap5  50wTrap6  50wTrap7  50wTrap8  50wTrap9 
    # 1826      2245       558       407      3174      1322 

set.seed(2565)
g10 <- buildSNNGraph(sce, k=10, use.dimred = 'HARMONY')
clust10 <- igraph::cluster_walktrap(g10)$membership
colData(sce)$k_10_Erik <-  paste0("10wTrap", factor(clust10))
table(colData(sce)$k_10_Erik)
    # 10wTrap1 10wTrap10 10wTrap11 10wTrap12 10wTrap13 10wTrap14 10wTrap15 10wTrap16 
    # 323        38      3218      1572      4254       216       149       499 
    # 10wTrap17 10wTrap18 10wTrap19  10wTrap2 10wTrap20 10wTrap21 10wTrap22 10wTrap23 
    # 202       169        97      1063       153       321        48       209 
    # 10wTrap24 10wTrap25 10wTrap26 10wTrap27 10wTrap28 10wTrap29  10wTrap3 10wTrap30 
    # 280        84       304        63        65        98      1326        65 
    # 10wTrap31 10wTrap32 10wTrap33 10wTrap34 10wTrap35 10wTrap36  10wTrap4  10wTrap5 
    # 22        37        39        17        38        17       138       179 
    # 10wTrap6  10wTrap7  10wTrap8  10wTrap9 
    # 686       337       371       385 

set.seed(14)
g5 <- buildSNNGraph(sce, k=5, use.dimred = 'HARMONY')
clust5 <- igraph::cluster_walktrap(g5)$membership
colData(sce)$k_5_Erik <- paste0("5wTrap", factor(clust5))
table(colData(sce)$k_5_Erik)
    # 5wTrap1  5wTrap10 5wTrap100 5wTrap101 5wTrap102 5wTrap103 5wTrap104 5wTrap105 
    # 61       275        20        10        14        23        14        11 
    # 5wTrap106 5wTrap107 5wTrap108  5wTrap11  5wTrap12  5wTrap13  5wTrap14  5wTrap15 
    # 7         7         8       337       413        97        70        48 
    # 5wTrap16  5wTrap17  5wTrap18  5wTrap19   5wTrap2  5wTrap20  5wTrap21  5wTrap22 
    # 94       142       127       133        47       112       120        33 
    # 5wTrap23  5wTrap24  5wTrap25  5wTrap26  5wTrap27  5wTrap28  5wTrap29   5wTrap3 
    # 120      1691       659        25       949        43        91        76 
    # 5wTrap30  ...

##### Approach 3: Mini-batch k Means (Steph Hicks Example) #####################
set.seed(789)

k_list <- seq(5, 20)

km_res <- lapply(k_list, function(k) {
  mbkmeans(sce, clusters = k, 
           batch_size = 500,
           reduceMethod = "HARMONY",
           calc_wcss = TRUE)
})

wcss <- sapply(km_res, function(x) sum(x$WCSS_per_cluster))

pdf(file = here("plots", "04_snRNA-seq", "06_Clustering", "mini_batch_klist_vs_wcss.pdf"),
    width = 5, height = 4)
plot(k_list, wcss, type = "b")
dev.off()
# interesting plot, going to save all k-means tho

# saving 
sce$kmeans5 <- paste0("mbk", km_res[[which(k_list==5)]]$Clusters)
sce$kmeans6 <- paste0("mbk", km_res[[which(k_list==6)]]$Clusters)
sce$kmeans7 <- paste0("mbk", km_res[[which(k_list==7)]]$Clusters)
sce$kmeans8 <- paste0("mbk", km_res[[which(k_list==8)]]$Clusters)
sce$kmeans9 <- paste0("mbk", km_res[[which(k_list==9)]]$Clusters)
sce$kmeans10 <- paste0("mbk", km_res[[which(k_list==10)]]$Clusters)
sce$kmeans11 <- paste0("mbk", km_res[[which(k_list==11)]]$Clusters)
sce$kmeans12 <- paste0("mbk", km_res[[which(k_list==12)]]$Clusters)
sce$kmeans13 <- paste0("mbk", km_res[[which(k_list==13)]]$Clusters)
sce$kmeans14 <- paste0("mbk", km_res[[which(k_list==14)]]$Clusters)
sce$kmeans15 <- paste0("mbk", km_res[[which(k_list==15)]]$Clusters)
sce$kmeans16 <- paste0("mbk", km_res[[which(k_list==16)]]$Clusters)
sce$kmeans17 <- paste0("mbk", km_res[[which(k_list==17)]]$Clusters)
sce$kmeans18 <- paste0("mbk", km_res[[which(k_list==18)]]$Clusters)
sce$kmeans19 <- paste0("mbk", km_res[[which(k_list==19)]]$Clusters)
sce$kmeans20 <- paste0("mbk", km_res[[which(k_list==20)]]$Clusters)

############# Pre-PLOTTING #####################################################
## NOTE: 
## 1) Getting rid of walkTrap 5 because we do not need 100+ grups 

# Fixing NAs in ct_Erik
sce$ct_Erik <- as.character(sce$ct_Erik)
sce$ct_Erik[is.na(sce$ct_Erik)] <- "Bukola"

## Fixing factors for sce$ct_Erik
sce$ct_Erik <- factor(sce$ct_Erik2, levels = c("Astro", "Endo", "Micro", "Oligo.1", "Oligo.2", "OPC.1", 
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
colorbyGroup <- c("louvain", "k_20_Erik", "k_50_Erik", "k_10_Erik", 
                  colnames(colData(sce)[grepl("kmeans", colnames(colData(sce)))]))
## Getting rid of walkTrap 5 because we do not need 100+ grups 
# [1] "louvain"   "k_20_Erik" "k_50_Erik" "k_10_Erik"  "kmeans5"  
# [7] "kmeans6"   "kmeans7"   "kmeans8"   "kmeans9"   "kmeans10"  "kmeans11" 
# [13] "kmeans12"  "kmeans13"  "kmeans14"  "kmeans15"  "kmeans16"  "kmeans17" 
# [19] "kmeans18"  "kmeans19"  "kmeans20"

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
         "Heatmap_Clusters_unfinalized4.pdf"), height = 7, width = 7)

additiongColorGroup <- c("k_20_Erik", "k_50_Erik", "k_10_Erik")

lapply(colorbyGroup, function(n) {
  # adding Rand index (index of similary on scale over 1)
  Rand <- signif(pairwiseRand(sce$ct_Erik, colData(sce)[, n], mode = "index"), 3)
  plot_cap <- paste("Rand Index Against Erik Cell Types:", Rand)
  
  linker <- linkClustersMatrix(sce$ct_Erik, colData(sce)[, n])
  
  # Not mark 0s as NAs because it messes with clustering.
  # linker[linker == 0] <- NA
  
  # heatmap 3  
  # Heatmap(linker,
  #         column_title = plot_cap,
  #         col = c("black", viridisLite::plasma(101)),
  #         name = "Corr"
  # )
  
  # heatmap 4
   pheatmap(linker,
         main = plot_cap,
   )
  
  
}) 

lapply(additiongColorGroup, function(n) {
  # adding Rand index (index of similary on scale over 1)
  Rand <- signif(pairwiseRand(sce$louvain, colData(sce)[, n], mode = "index"), 3)
  plot_cap <- paste("Rand Index Against Louvain Clusters:", Rand)
  
  linker <- linkClustersMatrix(sce$louvain, colData(sce)[, n])
  
  # Mark 0s as NAs
  linker[linker == 0] <- NA

  Heatmap(linker,
          column_title = plot_cap,
          col = viridisLite::plasma(101),
          na_col = "black",
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

## facet wrap by high Rand clusters and color with doublet score
## also change coloring scheme for clusters to show differences better



### SAVING OBJECT(S) ###########################################################
save(sce, file = here("processed-data", "04_snRNA-seq", "sce_objects", 
                      "sce_mid_clustering.Rdata"))


# bug check
table(sce$ct_Erik)

# Astro         Endo        Micro      Oligo.1      Oligo.2        OPC.1 
# 495          542          887          302           72          487 
# OPC.2        LHb.1        LHb.2        LHb.3        LHb.4        LHb.5 
# 197          317            5         1680          426           78 
# LHb.6        MHb.1        MHb.2 Neuron.Ambig  Thal.GABA.1  Thal.GABA.2 
# 220         1344          744          614          494          745 
# Thal.GABA.3      Thal.MD      Thal.PF     Thal.PVT       Bukola 
# 2632         4269           80          232          220 

# loading Erik's orignal cell types
load(here("processed-data", "08_snRNA-seq_Erik", "s3e_hb.rda"))


# Grabbing BrainNumber and Cell Type by Erik for each droplet
annoData <- data.frame(row.names = colnames(s3e.hb), "SampleID" = 
                         s3e.hb$sample_name, "ClusterID" = s3e.hb$cellType)
# Clearly identifying that each droplet is a nucleus
annoData$Barcode <- rownames(annoData)

# changing rownames
rownames(annoData) <- paste0(annoData$SampleID, "_", annoData$Barcode)

# dimensions for future reference [data from Erik]
dim(annoData)
# [1] 17529     3

# Creating a cell-type indicator in sce object's colData 
sce$ct_Erik2 <- NA
sce$ct_Erik2 <- annoData[colnames(sce),]$ClusterID

# checking
table(sce$ct_Erik2)

  # Astro         Endo        Micro      Oligo.1      Oligo.2        OPC.1 
  # 492           78          315         1674          425          220 
  # OPC.2 Neuron.Ambig        LHb.1        LHb.2        LHb.3        LHb.4 
  # 743            5         1341          743          613          494 
  # LHb.5        LHb.6        MHb.1        MHb.2  Thal.GABA.1  Thal.GABA.2 
  # 302           72          486          197         2626         4261 
  # Thal.GABA.3      Thal.MD      Thal.PF     Thal.PVT 
  # 79          231          220          541


## setting ct_Erik2 as ct_Erik (fixed in lines)

# checker: table(sce$ct_Erik, sce$ct_Erik2)
# getting rid of our secondary ct_Erik2 after using for fix in lines 158-163
# sce$ct_Erik2 <- NULL