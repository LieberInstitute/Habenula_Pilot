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

############### DOUBLET SCORE DROP #############################################
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


######## Table summary for drops by ct_Erik
ctErik_drop <- as.data.frame(table(sce$discard, sce$ct_Erik, by = sce$Sample))
names(ctErik_drop) <- c("T_F", "CellType", "BrNum", "Frequency")
ctErik_drop$BrNum <- as.character(ctErik_drop$BrNum)

barPlotErik <- function(BrNumber){
  plotter <- ctErik_drop[ctErik_drop$BrNum == BrNumber,]
  
  plotted <- ggplot(plotter, aes(x = CellType, y = Frequency, fill = T_F)) +
  geom_col(position = position_dodge()) +
  theme_bw() + 
  scale_fill_brewer(palette = "Greens") +
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
  lapply(unique(ctErik_drop$BrNum), barPlotErik)
dev.off()

############### Plotting before drops. #########################################
# per metric 
plotDropbySamp <- melt(toDropbySamp, id = "Sample")
plotDropbySamp$ToF <- ss(as.character(plotDropbySamp$variable), "_", 1)
plotDropbySamp$Metric <- sub("*._", "", as.character(plotDropbySamp$variable)) 

forSummaryTable <- function(BrNumber){
  plotting <- plotDropbySamp[plotDropbySamp$Sample == BrNumber,]
  plotting[plotting == "T"] <- "True"
  plotting[plotting == "F"] <- "False"
  plotting[plotting == "High_Mito"] <- "High Mito?"
  plotting[plotting == "Low_Lib"] <- "Low Lib?"
  plotting[plotting == "Low_Det_Feat"] <- "Low Det Fea?"
  
  plot <- ggplot(plotting,
       aes(x = Metric,
           y = value, 
           fill = ToF)) + 
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(title = BrNumber, fill = "Drop?", x = "Metric", y = element_blank()) +
    geom_text(aes(label = value), vjust  = -0.3, position = position_dodge(0.9), 
              size = 4, fontface = "bold") +
    scale_fill_brewer(palette = "Reds") +
    theme_minimal() + 
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 12, vjust = 0.2),
      axis.title.y = element_text(size = 12, vjust = 0.2),
      axis.text = element_text(face = "bold")
    )
  
  return(plot)
}

plottingDropbySamp <- lapply(unique(plotDropbySamp$Sample), forSummaryTable)
plottingDropbySamp[[8]] <- get_legend(plottingDropbySamp[[1]]) 

for (i in 1:7) {
  plottingDropbySamp[[i]] <- plottingDropbySamp[[i]] + theme(legend.position="none")
}

pdf(file = here("plots", "04_snRNA-seq", "03_qc_metter_plots", "bar_plots_tf_drops.pdf"), 
    width = 9, height = 16)
      do.call("grid.arrange", c(plottingDropbySamp, ncol = 2))
dev.off()

# Overall drop info bar plot summary 
overallDropDF$Sample <- rownames(overallDropDF) 
overallDropDF <- melt(overallDropDF, by = Sample)
overallDropDF$variable <- recode(overallDropDF$variable,  "totalKeep" = 
                                 "Total Post-Drop", "totalDiscard" = "Dropped",
                                 "totalDroplets" = "Total Pre-Drop")
overallDropDF$variable <- factor(overallDropDF$variable, levels = 
                                   c("Total Pre-Drop", "Total Post-Drop", "Dropped"))

plotOverall <- ggplot(overallDropDF,
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

pdf(file = here("plots", "04_snRNA-seq", "j03_qc_metter_plots", "overall_drop_barplot.pdf"), 
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
    ggtitle("Detected Features")
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


## Do doublet scores before dropping doublet scores and then make a table regarding
# how many cells are dropped after QC. - Erik 
