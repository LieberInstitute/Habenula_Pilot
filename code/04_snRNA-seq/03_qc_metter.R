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

############### Plotting before drops. #########################################
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

# Plotting before and after dropping info
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
    geom_bar(stat="identity", position = position_dodge()) +
    labs(title = BrNumber, fill = "Drop?", x = "Metric", y = element_blank()) +
    geom_text(aes(label = value), vjust  = -0.3, position = position_dodge(0.9), 
              size = 4, fontface = "bold") +
    scale_fill_brewer(palette = "Reds") +
    theme_minimal() + 
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
    )
  
  return(plot)
}

plottingDropbySamp <- lapply(unique(plotDropbySamp$Sample), forSummaryTable)
plottingDropbySamp[[8]] <- get_legend(plottingDropbySamp[[1]])

pdf(file = here("plots", "04_snRNA-seq", "03_qc_metter_plots", "bar_plots_tf_drops.pdf"))
  do.call("grid.arrange", c(plottingDropbySamp, ncol = 2))
dev.off()


# recording sce pre drop
dim_predrop = dim(sce)
dim_predrop
    # 36601 19802

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


############### DROPPING TRUES #################################################
# Dropping high mito, low library size, and low features
sce <- sce[, sce$discard == FALSE]

# Dropping 0 count genes across samples
sce <- sce[rowSums(counts(sce)) > 0, ]

# recording sce post drop
dim_postmitodrop = dim(sce)
dim_postmitodrop
# 33848 17082

pdf(file = here("plots", "04_snRNA-seq", "03_qc_metter_plots", "dist_vplot_post_drop.pdf"))
  plotColData(sce, x = "Sample", y = "subsets_Mito_percent") + ggtitle("Mito Precent")
  
  plotColData(sce, x = "Sample", y = "sum") + ggtitle("Library Size")
  
  plotColData(sce, x = "Sample", y = "detected") + ggtitle("Detected Features")
dev.off()

sce_hb_postQC <- sce 
rm(sce)

save(sce_hb_postQC, file = here("processed-data", "04_snRNA-seq", "sce_objects", 
                                "sce_hb_postQC.Rdata"))

############### CLUSTER ANNOTATIOS {Erik's} ####################################
load(here("processed-data", "08_snRNA-seq_Erik", "s3e_hb.rda")) 

annoData <- data.frame(row.names = colnames(s3e.hb), "SampleID" = 
                         s3e.hb$sample_name, "ClusterID" = s3e.hb$cellType)




# Questions for Louise:
# 1) Should we run qc per sample before combining them for the sce object

## 1a: No, add column for batch and then add argument to qc mito 

# 2) Why does the sce object look like that. Sample and Barcode only in colData vs 
# that of Erik's s3e.hb.
## 2a: Change the Sample to the Br name it is asscociated with 

# 3) should we use josh's or erik's way of grabbing mito gennes

# 4) aren't we supposed to run addpercellQC because I can't find sum, detected, 
# etc for plotting
# 4a: that was unfilitered. just add all info to same sce object and 

# LATER: 3) Drop any genes with 0 counts in any nuclei 
