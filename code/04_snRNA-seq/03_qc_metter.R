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

# Loading filtered pre-QC sce object for all 7 Habenula samples
load(here("processed-data", "04_snRNA-seq", "sce_objects", "sce_hb_preQC.Rdata"))
sce <- sce_hb_preQC
rm(sce_hb_preQC)

# QC approach based on Workflow 3.3 in OSCA:
# http://bioconductor.org/books/3.16/OSCA.workflows/unfiltered-human-pbmcs-10x-genomics.html#quality-control-2

############### HIGH MITO DROP 
# --------------- Way 1:
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

# Identifying high mito reads
sce$high_mito <- isOutlier(sce$subsets_Mito_percent, nmads=3, type="higher", 
                            batch = sce$Sample)
table(sce$high_mito)
    # high_mito
    # FALSE  TRUE 
    # 17702  2100 

# recording sce pre drop
dim_premitodrop = dim(sce)
dim_premitodrop
    # 36601 19802

# Dropping high mito
sce <- sce[, sce$high_mito == FALSE]

# recording sce post drop
dim_postmitodrop = dim(sce)
dim_postmitodrop
    # 36601 17245

# for plotting
sce$discard = sce$high_mito

# Plotting qc metric distribution
pdf(file = here("plots", "04_snRNA-seq", "03_qc_metter_plots", "high_mito_dist.R"))
  plotColData(sce, x = "Sample", y = "subsets_Mito_percent", colour_by = "high_mito") +
    ggtitle("Mito Precent")
dev.off()

############### LOW LIBRARY SIZE DROP



############### LOW DETECTED FEATURES DROP



# Notes (from Erik):
# 1) Compare (location=="MT") vs grep("^MT-")
# 2) Drop based on qc metrics
# 3) Drop any genes with 0 counts in any nuclei 

# Getting cluster annotations from Erik's sce object
load(here("processed-data", "08_snRNA-seq_Erik", "s3e_hb.rda")) 

annoData <- data.frame(row.names = colnames(s3e.hb), "SampleID" = s3e.hb$sample_name,
                         "ClusterID" = s3e.hb$cellType)




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