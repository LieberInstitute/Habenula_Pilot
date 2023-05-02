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
pdf(here(plot_dir, "sce_pre_QC_violin_plots_by_Run.pdf"))

plotColData(sce, x = "Sample", y = "subsets_Mito_percent", colour_by = "high_mito") + 
  ggtitle("Mito Precent")

plotColData(sce, x = "Sample", y = "sum", colour_by = "lowLib") +
  ggtitle("Library Size")

plotColData(sce, x = "Sample", y = "detected", colour_by = "lowDetecFea") +
  ggtitle("Detected Features")

dev.off()




# .



print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()