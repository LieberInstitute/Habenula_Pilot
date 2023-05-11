## May 2, 2023 - Bukola Ajanaku
# sce data pre-qc on violin plots to show threshold for drops. 
# qrsh -l mem_free=30G,h_vmem=30G

# loading relevant libraries
library("SingleCellExperiment")
library("here")
library("ggplot2")
library("EnsDb.Hsapiens.v86")
library("scuttle")
library("jaffelab")
library("sessioninfo")
library("scater")
library("PupillometryR")

# loading pre-QC sce object
load(here("processed-data", "99_paper_figs",  "sce_objects", "pre_QC_sce.Rdata"),
     verbose = TRUE)
# sce_hb_preQC
sce <- sce_hb_preQC
rm(sce_hb_preQC)

# creating plot directory
plot_dir <- here("plots", "99_paper_figs", "04_QC_Violin.R")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

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


################ PLOTTING PER METRIC ###########################################
pd <- as.data.frame(colData(sce))
sce$lowLib <- as.factor(sce$lowLib)
sce$Sample <- as.factor(sce$Sample)

## high_mito
pdf(here(plot_dir, "official_Violin_QC_High_Mito.pdf"), width = 8)
  ggplot(pd, aes(x = Sample, y = subsets_Mito_percent)) +
    geom_jitter(aes(color = high_mito), position = 
                 position_jitter(seed = 1, width = 0.01)) +
    labs(y = "Mitochodrial Percent", colour = "High Mito?") +
    geom_violin(alpha = 0.4, position = "identity", scale = 3) +
    scale_y_log10() +
    geom_violin(
                # data = pd[pd$high_mito == FALSE, ],
                alpha = 0.4, position = "identity",
                scale = 3, aes(fill = high_mito)) + 
                guides(fill = FALSE) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          axis.title.x = element_blank()) +
    theme_bw() + 
    scale_y_log10()
dev.off()

# low libray size
pdf(here(plot_dir, "official_Violin_QC_Low_Lib.pdf"), width = 8)
  ggplot(pd, aes(x = Sample, y = sum, group = Sample, fill = lowLib)) +
    geom_point(aes(color = lowLib, group = Sample)) +
    geom_violin(data = pd, alpha = 0.4, position = "identity", 
                inherit.aes = TRUE) +
    labs(y = "Library Size", colour = "Low Lib Size?") +
    scale_y_log10()

dev.off()

pdf(here(plot_dir, "test.pdf"))
  plotColData(sce, x = "Sample", y="sum", colour_by="lowLib") +
    scale_y_log10() + ggtitle("Total count") + aes(group = pd$Sample)
dev.off()



print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()