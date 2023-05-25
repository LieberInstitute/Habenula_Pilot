# January 12, 2023 - Bukola Ajanaku
# Applying QC metrics to filtered sce objects of the 7 Habenula samples.
# qrsh -l mem_free=50G,h_vmem=50G

# Loading relevant libraries
library("SingleCellExperiment")
library("jaffelab")
library("VariantAnnotation")
library("here")
library("ggplot2")
library("ggrepel")
library("scater")
library("batchelor")
library("scran")
library("scry")
library("uwot")
library("DropletUtils")
library("Rtsne")
library("gridExtra")
library("EnsDb.Hsapiens.v86")
library("reshape")
library("cowplot")
library("dplyr")
library("scDblFinder")
library("sessioninfo")

# Loading filtered pre-QC sce object for all 7 Habenula samples
load(here("processed-data", "04_snRNA-seq", "sce_objects", "sce_hb_preQC.Rdata"))
sce <- sce_hb_preQC
rm(sce_hb_preQC)

# QC approach based on Workflow 3.3 in OSCA:
# http://bioconductor.org/books/3.16/OSCA.workflows/unfiltered-human-pbmcs-10x-genomics.html#quality-control-2

############### FINDING HIGH MITO ##############################################
# MADs approach for high mito droplets (indicates nuclei where cytoplasm wasn't
# successfully stripped):
set.seed(777)

# Defining locations
location <- mapIds(EnsDb.Hsapiens.v86, keys = rowData(sce)$ID, 
                   column = "SEQNAME", keytype = "GENEID")
# Unable to map 2611 of 36601 requested IDs.

# Identifying mito reads
stats <- perCellQCMetrics(sce, subsets = list(Mito = which(location=="MT")))

# Changing Sample ID names from locations 
sce$path <- sce$Sample
sce$Sample <- ss(sce$Sample, "/", 9)

# Binding stats to colData of sce
colData(sce) <- cbind(colData(sce), stats)

############### HIGH MITO ######################################################
sce$high_mito <- isOutlier(sce$subsets_Mito_percent, nmads=3, type="higher", 
                            batch = sce$Sample)
table(sce$high_mito, by = sce$Sample)
    #       Br1092 Br1204 Br1469 Br1735 Br5555 Br5558 Br5639
    # FALSE   3352   1391   2185   3199   3600    801   3174
    # TRUE     270    274    314    544    305    135    258

############### LOW LIBRARY SIZE ###############################################
sce$lowLib <- isOutlier(sce$sum, log=TRUE, type="lower", batch = sce$Sample)
table(sce$lowLib, by = sce$Sample)
    #       Br1092 Br1204 Br1469 Br1735 Br5555 Br5558 Br5639
    # FALSE   3430   1445   2414   3633   3772    847   3348
    # TRUE     192    220     85    110    133     89     84

############### LOW DETECTED FEATURES ##########################################
sce$lowDetecFea <- isOutlier(sce$detected, log=TRUE, type="lower", batch = sce$Sample)
table(sce$lowDetecFea, by = sce$Sample)
    #       Br1092 Br1204 Br1469 Br1735 Br5555 Br5558 Br5639
    # FALSE   3394   1334   2347   3437   3766    826   3289
    # TRUE     228    331    152    306    139    110    143

############### DOUBLET SCORE ##################################################
set.seed(234)

colData(sce)$doubletScore <- NA

for (i in splitit(sce$Sample)) {
  sce_temp <- sce[, i]
  ## To speed up, run on sample-level top-HVGs - just take top 1000
  normd <- logNormCounts(sce_temp)
  geneVar <- modelGeneVar(normd)
  topHVGs <- getTopHVGs(geneVar, n = 1000)
  
  dbl_dens <- computeDoubletDensity(normd, subset.row = topHVGs)
  colData(sce)$doubletScore[i] <- dbl_dens
}

summary(sce$doubletScore)
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    # 0.0000  0.1484  0.3449  0.6870  0.7550 31.9605 

### Plotting (gotten from:  
# https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/code/03_build_sce/03_droplet_qc.R)

# plotting per Sample and indicating line of doublet score 5. Any higher is not good.
dbl_df <- colData(sce) %>%
  as.data.frame() %>%
  select(Sample, doubletScore)

dbl_box_plot <- dbl_df %>%
  ggplot(aes(x = reorder(Sample, doubletScore, FUN = median), y = doubletScore)) +
  geom_boxplot() +
  labs(x = "Sample") +
  geom_hline(yintercept = 5, color = "red", linetype = "dashed") +
  coord_flip() +
  theme_bw()

png(filename = here("plots", "04_snRNA-seq", "03_qc_metter_plots",
                    "doublet_scores_boxplot.png"))
  dbl_box_plot
dev.off()

# plotting density for doublet score 
dbl_density_plot <- dbl_df %>%
  ggplot(aes(x = doubletScore)) +
  geom_density() +
  labs(x = "doublet score") +
  facet_grid(Sample ~ .) +
  theme_bw()

png(filename = here("plots", "04_snRNA-seq", "03_qc_metter_plots",
                    "doublet_scores_density.png") , width = 700, height = 600)
  dbl_density_plot
dev.off()

table(sce$doubletScore >= 5)
    # FALSE  TRUE 
    # 19585   217 

# returning table regarding drop by doublet score
dbl_df %>%
  group_by(Sample) %>%
  summarize(
    median = median(doubletScore),
    q95 = quantile(doubletScore, .95),
    drop = sum(doubletScore >= 5),
    drop_precent = 100 * drop / n()
  )
    # Sample median   q95  drop drop_precent
    # <chr>   <dbl> <dbl> <int>        <dbl>
    #   1 Br1092  0.355  1.56    35       0.966 
    # 2 Br1204  0.556  1.32     1       0.0601
    # 3 Br1469  0.290  1.80    27       1.08  
    # 4 Br1735  0.299  1.76    45       1.20  
    # 5 Br5555  0.226  1.11    77       1.97  
    # 6 Br5558  0.958  2.52     0       0     
    # 7 Br5639  0.494  2.42    32       0.932 

############### LOADING ERIK'S ANNOTATED CLUSTERS ##############################
load(here("processed-data", "08_snRNA-seq_Erik", "s3e_hb.rda"))


# CHECKING ERIKS FOR TEST
test <- "TCGCTCACAAATAGCA-1"
table(duplicated(colnames(s3e.hb))) #none 
colData(sce)[sce$Barcode == test,] # Erik has test for 1204 and not 1092
colData(sce)[sce$Barcode == test,]$path


# Grabbing BrainNumber and Cell Type by Erik for each droplet
annoData <- data.frame(row.names = colnames(s3e.hb), "SampleID" = 
                         s3e.hb$sample_name, "ClusterID" = s3e.hb$cellType)
# Clearly identifying that each droplet is a nucleus
annoData$Barcode <- rownames(annoData)

# dimensions for future reference [data from Erik]
dim(annoData)
    # 17529     3

# dimesions for future reference [undropped sce, no qc applied just yet]
dim(colData(sce))  
  #  36601 19802

# Creating a cell-type indicator in sce object's colData 
sce$ct_Erik <- NA
sce$ct_Erik <- annoData$ClusterID[match(sce$Barcode, annoData$Barcode)]

# Sanity check!
table(is.na(sce$ct_Erik))
    # FALSE  TRUE 
    # 17000  2802

## Updating annoData with my sce information
# Barcode only in Barcode column
rownames(annoData) <- NULL

# Defaulting all barcodes to No and only changing if confirmed to be in my sce object
annoData$inBukolas <- "No"

# Confirming if present in my sce object
annoData[annoData$Barcode %in% sce$Barcode,]$inBukolas <- "Yes"

# grabbing barcode data that are not present in Erik's
notinEriks <- colData(sce[, !(sce$Barcode %in% annoData$Barcode)])

# subsetting to only relevant data per barcode
notinEriks <- notinEriks[, c("Sample", "Barcode")]
rownames(notinEriks) <- NULL
colnames(notinEriks)[1] <- "SampleID"

# Adding columns to notinEriks pre-rbind
notinEriks$inBukolas <- "Yes"
notinEriks$ClusterID <- "New.TBD"
notinEriks <- as.data.frame(notinEriks)
# continuing with Erik df after reorganization of drop info.

############### Re-organizing drop info. #######################################
# recording sce pre drop
dim_predrop = dim(sce)
dim_predrop
   # 36601 19802

# total we want to drop 
sce$discard <- sce$high_mito | sce$lowLib | sce$lowDetecFea
table(sce$discard)
    # FALSE  TRUE 
    # 17082  2720 

# data frame of qc met stats 
# Per Metric
toDropbySamp <- t(table(sce$high_mito, by = sce$Sample))
colnames(toDropbySamp) <- c("F_High_Mito", "T_High_Mito") 

toDropbySamp <- cbind(toDropbySamp, t(table(sce$lowLib, by = sce$Sample)))
colnames(toDropbySamp)[3:4] <- c("F_Low_Lib", "T_Low_Lib") 

toDropbySamp <- cbind(toDropbySamp, t(table(sce$lowDetecFea, by = sce$Sample)))
colnames(toDropbySamp)[5:6] <- c("F_Low_Det_Feat", "T_Low_Det_Feat") 

toDropbySamp <- as.data.frame(toDropbySamp)
toDropbySamp$Sample <- rownames(toDropbySamp)

# Overall 
overallDropDF <- t(table(sce$discard, by = sce$Sample))
colnames(overallDropDF) <- c("totalKeep", "totalDiscard") 
overallDropDF <- cbind(overallDropDF, rowSums(overallDropDF))
colnames(overallDropDF)[3] <- "totalDroplets"
overallDropDF <- as.data.frame(overallDropDF)

######## Continuing ERIK'S ANNOTATED CLUSTERS DF ###############################
# Current bar plots subsetting sce by annoData and adding new nuc from annoData:
      # # Combining my kept barcodes to Erik's df of annotated barcordes
      # annoData <- rbind(annoData, notinEriks)
      # table(annoData$inBukolas)
      # # No   Yes 
      # # 570 19761

# Goal: saving df regarding sce object from Erik, indicating samples he kept that I didn't
# and my samples that made the cut for my sce object that aren't in his.

# grabbing from my sce
combinedData <- colData(sce)[c("Sample", "ct_Erik", "high_mito", "lowDetecFea", "lowLib", "discard")]
colnames(combinedData)[1] <- "SampleID"
colnames(combinedData)[2] <- "ClusterID"
combinedData$Barcode <- rownames(combinedData)
rownames(combinedData) <- NULL
# make sure to prep for combination
combinedData$inBukolas <- "Yes"
combinedData <- as.data.frame(combinedData)

# Subsetting annoData for only nuc that are not already in my sce object
ErikOnly <- annoData[!(annoData$Barcode %in% combinedData$Barcode),]

# prepping annoData for combination
ErikOnly$high_mito <- "ErikKept"
ErikOnly$lowDetecFea <- "ErikKept"
ErikOnly$lowLib <- "ErikKept"
ErikOnly$discard <- "ErikKept"
ErikOnly$inBukolas <- "No" 
  # Should be 570 Nos

# Sanity check
dim(ErikOnly)
  # 570   8
dim(combinedData) #pre-combination
  # 19802     8

# officially creating combinedData
combinedData <- rbind(combinedData, ErikOnly)
dim(combinedData)
  # 20372     8

# Adding the info for a possible new cluster made by me
forNewCluster <- annoData[annoData$ClusterID == "New.TBD",]$Barcode
length(forNewCluster)
  # 2802

### OFFICIALLY COMPLETE COMBINED DATA (annoData including Bukola new samples) ***
combinedData[combinedData$Barcode %in% forNewCluster,]$ClusterID <- "New.TBD"

table(combinedData$inBukolas)
    # No   Yes 
    # 570 19802 

table(combinedData$ClusterID)
    # Astro         Endo        Micro      Oligo.1      Oligo.2        OPC.1 
    # 612          118          381         1807          427          278 
    # OPC.2 Neuron.Ambig        LHb.1        LHb.2        LHb.3        LHb.4 
    # 767           79         1361          767          619          516 
    # LHb.5        LHb.6        MHb.1        MHb.2  Thal.GABA.1  Thal.GABA.2 
    # 302           74          497          197         2916         4709 
    # Thal.GABA.3      Thal.MD      Thal.PF     Thal.PVT      New.TBD 
    # 95          262          239          547         2802

############################ Plotting ##########################################
# ctErik_drop <- table(sce$discard, sce$ct_Erik, by = sce$Sample)
ctErik_drop <- table(combinedData$discard, combinedData$ClusterID, by = combinedData$SampleID)
ctErik_drop <- as.data.frame(ctErik_drop)
  # True = Drop for me 
  # False = Do not drop for me
  # Erik Kept = What Erik Kept but I didn't
names(ctErik_drop) <- c("NotKept", "ClusterID", "SampleID", "Frequency")

barPlotErik <- function(BrNumber){
  plotter <- ctErik_drop[ctErik_drop$SampleID == BrNumber,]
  
  plotted <- ggplot(plotter, 
                    aes(x = ClusterID, y = Frequency, fill = NotKept)) +
  geom_col(position = position_dodge()) +
  theme_bw() + 
  scale_fill_brewer(palette="Dark2") +
  labs(title = BrNumber, fill = "Drop?", x = "Cell Type", y = element_blank()) +
  geom_text(
    aes(label = Frequency),
    colour = "darkblue", 
    size = 3 , 
    position = position_dodge(.9)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(0, 250))
  
  return(plotted)
}


pdf(file = here("plots", "04_snRNA-seq", "03_qc_metter_plots", "sce_qc_cellType_Erik.pdf"),
    width = 14, height = 7)
  lapply(unique(ctErik_drop$SampleID), barPlotErik)
dev.off()

############### TAKE 2: Re-organizing drop info. ###############################
# Using combinedData df 

# data frame of qc met stats 
# Per Metric *Need to split Erik Kept by type*
newDropbySamp <- t(table(combinedData$high_mito, by = combinedData$SampleID))
colnames(newDropbySamp) <- c("ErikKept", "F_High_Mito", "T_High_Mito") 

newDropbySamp <- cbind(newDropbySamp, t(table(combinedData$lowLib, by = combinedData$SampleID)))
colnames(newDropbySamp)[4:6] <- c("ErikKept", "F_Low_Lib", "T_Low_Lib") 

newDropbySamp <- cbind(newDropbySamp, t(table(combinedData$lowDetecFea, by = combinedData$SampleID)))
colnames(newDropbySamp)[7:9] <- c("ErikKept", "F_Low_Det_Feat", "T_Low_Det_Feat") 

newDropbySamp <- as.data.frame(newDropbySamp)
newDropbySamp$Sample <- rownames(newDropbySamp)

# Overall 
newOverallDropDF <- t(table(combinedData$discard, by = combinedData$Sample))
colnames(newOverallDropDF) <- c("ErikKept","totalKeep", "totalDiscard") 
newOverallDropDF <- cbind(newOverallDropDF, rowSums(newOverallDropDF))
colnames(newOverallDropDF)[4] <- "totalDroplets"
newOverallDropDF <- as.data.frame(newOverallDropDF)
newOverallDropDF$Sample <- rownames(newDropbySamp)

############### Plotting before drops. #########################################
# per metric 
newPlotDropbySamp <- newDropbySamp
newPlotDropbySamp <- newPlotDropbySamp[-c(4,7)]
rownames(newPlotDropbySamp) <- NULL
newPlotDropbySamp <- melt(newPlotDropbySamp, by = newPlotDropbySamp$Sample)
colnames(newPlotDropbySamp) <- c("Sample", "NotKept", "Value")

newPlotDropbySamp$ToF <- ss(as.character(newPlotDropbySamp$NotKept), "_", 1)
newPlotDropbySamp$Metric <- sub("*._", "", as.character(newPlotDropbySamp$NotKept)) 
newPlotDropbySamp$Sample <- as.character(newPlotDropbySamp$Sample)

forSummaryTable <- function(BrNumber){
  plotting <- newPlotDropbySamp[newPlotDropbySamp$Sample == BrNumber,]
  plotting[plotting == "T"] <- "True"
  plotting[plotting == "F"] <- "False"
  plotting[plotting == "High_Mito"] <- "High Mito?"
  plotting[plotting == "Low_Lib"] <- "Low Lib?"
  plotting[plotting == "Low_Det_Feat"] <- "Low Det Fea?"
  
  plot <- ggplot(plotting,
       aes(x = Metric,
           y = Value, 
           fill = ToF)) + 
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(title = BrNumber, fill = "Drop?", x = "Metric", y = element_blank()) +
    geom_text(aes(label = Value), vjust  = -0.3, position = position_dodge(0.9), 
              size = 4, fontface = "bold") +
    scale_fill_brewer(palette = "Accent") +
    theme_minimal() + 
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 12, vjust = 0.2),
      axis.title.y = element_text(size = 12, vjust = 0.2),
      axis.text = element_text(face = "bold")
    )
  
  return(plot)
}

plottingDropbySamp <- lapply(unique(newPlotDropbySamp$Sample), forSummaryTable)
plottingDropbySamp[[8]] <- get_legend(plottingDropbySamp[[1]]) 

for (i in 1:7) {
  plottingDropbySamp[[i]] <- plottingDropbySamp[[i]] + theme(legend.position="none")
}

pdf(file = here("plots", "04_snRNA-seq", "03_qc_metter_plots", "bar_plots_tf_drops_update.pdf"), 
    width = 9, height = 16)
      do.call("grid.arrange", c(plottingDropbySamp, ncol = 2))
dev.off()

# Overall drop info bar plot summary 
newOverallDropDF$Sample <- rownames(overallDropDF) 
newOverallDropDF <- melt(newOverallDropDF, by = Sample)
newOverallDropDF$variable <- recode(newOverallDropDF$variable,  "totalKeep" = 
                                 "Total Post-Drop", "totalDiscard" = "Dropped",
                                 "totalDroplets" = "Total Pre-Drop", 
                                 "ErikKept" = "Erik Kept It")
newOverallDropDF$variable <- factor(newOverallDropDF$variable, levels = 
                                   c("Total Pre-Drop", "Total Post-Drop", "Dropped",
                                     "Erik Kept It"))

plotOverall <- ggplot(newOverallDropDF,
               aes(x = Sample,
                   y = value, 
                   fill = variable)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Dropping Metrics", fill = "Sum Type", x = "Sample", y = element_blank()) +
  geom_text(aes(label = value), vjust  = -0.2, position = position_dodge(0.9), 
            size = 4, fontface = "bold") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12, vjust = 0.2),
    axis.title.y = element_text(size = 12, vjust = 0.2),
    axis.text = element_text(face = "bold")
  )

pdf(file = here("plots", "04_snRNA-seq", "03_qc_metter_plots", "overall_drop_barplot_update.pdf"), 
    width = 8, height = 6)
  plotOverall
dev.off()

# Plotting qc metric distribution
pdf(file = here("plots", "04_snRNA-seq", "03_qc_metter_plots", "dist_vplot_per_met.pdf"))
  plotColData(sce, x = "Sample", y = "subsets_Mito_percent", colour_by = "high_mito") + 
  ggtitle("Mito Precent")
  
  plotColData(sce, x = "Sample", y = "sum", colour_by = "lowLib") +
    ggtitle("Library Size")
  
  plotColData(sce, x = "Sample", y = "detected", colour_by = "lowDetecFea") +
    ggtitle("Detected Features")
dev.off()

pdf(file = here("plots", "04_snRNA-seq", "03_qc_metter_plots", "dist_vplot_by_all_discard.pdf"))
 plotColData(sce, x = "Sample", y = "subsets_Mito_percent", colour_by = "discard") +
    ggtitle("Mito Precent")
  
  plotColData(sce, x = "Sample", y = "sum", colour_by = "discard") +
    ggtitle("Library Siz")
  
  plotColData(sce, x = "Sample", y = "detected", colour_by = "discard") +
    ggtitle("Detected Features") + 
dev.off()


############### DROPPING TRUES & PLOTTING ######################################
# Dropping high mito, low library size, and low features
sce <- sce[, sce$discard == FALSE]

# Dropping 0 count genes across samples
sce <- sce[rowSums(counts(sce)) > 0, ]

# recording sce post drop
dim_postmitodrop = dim(sce)
dim_postmitodrop
    # 33848 17082

#### Plotting post drop!

pdf(file = here("plots", "04_snRNA-seq", "03_qc_metter_plots", "dist_vplot_post_drop.pdf"))
  plotColData(sce, x = "Sample", y = "subsets_Mito_percent") + ggtitle("Mito Precent")
  
  plotColData(sce, x = "Sample", y = "sum") + ggtitle("Library Size")
  
  plotColData(sce, x = "Sample", y = "detected") + ggtitle("Detected Features")
dev.off()

sce_hb_postQC <- sce 
rm(sce)



############### Saving sce objects made here. ##################################
save(sce_hb_postQC, file = here("processed-data", "04_snRNA-seq", "sce_objects", 
                                "sce_hb_postQC.Rdata"))

save(annoData, file = here("processed-data", "04_snRNA-seq", "sce_objects", 
                                "annoData.Rdata"))

# Saving for future reference
save(ctErik_drop, file = here("processed-data", "04_snRNA-seq", "sce_objects", 
                              "ctErik_drop.Rdata"))
# Saving combined info regard sce object from Erik, indicating sample he kept that I didn't
# and my samples that made the cut for my sce object that aren't in his.
save(combinedData, file = here("processed-data", "04_snRNA-seq", "sce_objects", 
                              "combinedData.Rdata"))

### sessionInfo()

# R version 4.2.3 Patched (2023-04-07 r84211)
# Platform: x86_64-conda-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)
# 
# Matrix products: default
# BLAS:   /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/lib/libRblas.so
# LAPACK: /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/lib/libRlapack.so
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices datasets  utils     methods  
# [8] base     
# 
# other attached packages:
#   [1] sessioninfo_1.2.2           scDblFinder_1.12.0         
# [3] dplyr_1.1.1                 cowplot_1.1.1              
# [5] reshape_0.8.9               EnsDb.Hsapiens.v86_2.99.0  
# [7] ensembldb_2.22.0            AnnotationFilter_1.22.0    
# [9] GenomicFeatures_1.50.4      AnnotationDbi_1.60.2       
# [11] gridExtra_2.3               Rtsne_0.16                 
# [13] DropletUtils_1.18.1         uwot_0.1.14                
# [15] Matrix_1.5-4                scry_1.10.0                
# [17] scran_1.26.2                batchelor_1.14.1           
# [19] scater_1.26.1               scuttle_1.8.4              
# [21] ggrepel_0.9.3               ggplot2_3.4.2              
# [23] here_1.0.1                  VariantAnnotation_1.44.1   
# [25] Rsamtools_2.14.0            Biostrings_2.66.0          
# [27] XVector_0.38.0              jaffelab_0.99.32           
# [29] rafalib_1.0.0               SingleCellExperiment_1.20.1
# [31] SummarizedExperiment_1.28.0 Biobase_2.58.0             
# [33] GenomicRanges_1.50.2        GenomeInfoDb_1.34.9        
# [35] IRanges_2.32.0              S4Vectors_0.36.2           
# [37] BiocGenerics_0.44.0         MatrixGenerics_1.10.0      
# [39] matrixStats_0.63.0          colorout_1.2-2             
# 
# loaded via a namespace (and not attached):
#   [1] BiocFileCache_2.6.1       plyr_1.8.8               
# [3] igraph_1.4.2              lazyeval_0.2.2           
# [5] splines_4.2.3             BiocParallel_1.32.6      
# [7] digest_0.6.31             viridis_0.6.2            
# [9] fansi_1.0.4               magrittr_2.0.3           
# [11] memoise_2.0.1             BSgenome_1.66.3          
# [13] ScaledMatrix_1.6.0        cluster_2.1.4            
# [15] limma_3.54.2              R.utils_2.12.2           
# [17] prettyunits_1.1.1         colorspace_2.1-0         
# [19] blob_1.2.4                rappdirs_0.3.3           
# [21] jsonlite_1.8.4            crayon_1.5.2             
# [23] RCurl_1.98-1.12           glue_1.6.2               
# [25] gtable_0.3.3              gargle_1.3.0             
# [27] zlibbioc_1.44.0           DelayedArray_0.24.0      
# [29] BiocSingular_1.14.0       Rhdf5lib_1.20.0          
# [31] HDF5Array_1.26.0          scales_1.2.1             
# [33] DBI_1.1.3                 edgeR_3.40.2             
# [35] Rcpp_1.0.10               viridisLite_0.4.1        
# [37] progress_1.2.2            dqrng_0.3.0              
# [39] bit_4.0.5                 rsvd_1.0.5               
# [41] ResidualMatrix_1.8.0      metapod_1.6.0            
# [43] httr_1.4.5                RColorBrewer_1.1-3       
# [45] pkgconfig_2.0.3           XML_3.99-0.14            
# [47] R.methodsS3_1.8.2         dbplyr_2.3.2             
# [49] locfit_1.5-9.7            utf8_1.2.3               
# [51] tidyselect_1.2.0          rlang_1.1.0              
# [53] munsell_0.5.0             tools_4.2.3              
# [55] cachem_1.0.7              xgboost_1.7.5.1          
# [57] cli_3.6.1                 generics_0.1.3           
# [59] RSQLite_2.3.1             stringr_1.5.0            
# [61] fastmap_1.1.1             yaml_2.3.7               
# [63] bit64_4.0.5               fs_1.6.1                 
# [65] purrr_1.0.1               KEGGREST_1.38.0          
# [67] nlme_3.1-162              sparseMatrixStats_1.10.0 
# [69] R.oo_1.25.0               xml2_1.3.3               
# [71] biomaRt_2.54.1            compiler_4.2.3           
# [73] beeswarm_0.4.0            filelock_1.0.2           
# [75] curl_5.0.0                png_0.1-8                
# [77] tibble_3.2.1              statmod_1.5.0            
# [79] stringi_1.7.12            lattice_0.20-45          
# [81] bluster_1.8.0             ProtGenerics_1.30.0      
# [83] vctrs_0.6.1               pillar_1.9.0             
# [85] lifecycle_1.0.3           rhdf5filters_1.10.1      
# [87] BiocNeighbors_1.16.0      data.table_1.14.8        
# [89] bitops_1.0-7              irlba_2.3.5.1            
# [91] rtracklayer_1.58.0        R6_2.5.1                 
# [93] BiocIO_1.8.0              vipor_0.4.5              
# [95] codetools_0.2-19          MASS_7.3-58.2            
# [97] rhdf5_2.42.1              rprojroot_2.0.3          
# [99] rjson_0.2.21              withr_2.5.0              
# [101] GenomicAlignments_1.34.1  GenomeInfoDbData_1.2.9   
# [103] parallel_4.2.3            hms_1.1.3                
# [105] grid_4.2.3                beachmat_2.14.2          
# [107] DelayedMatrixStats_1.20.0 segmented_1.6-2          
# [109] googledrive_2.1.0         ggbeeswarm_0.7.1         
# [111] restfulr_0.0.15          
 








