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

# Prepping sce object to be scanned for mito labels
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)

# Defining locations
location <- mapIds(EnsDb.Hsapiens.v86, keys = rowData(sce)$ID, 
                   column = "SEQNAME", keytype = "GENEID")
# Unable to map 2611 of 36601 requested IDs.

# Identifying mito reads
stats <- perCellQCMetrics(sce, subsets = list(Mito = which(location=="MT")))

# Identifying high mito reads
high_mito <- isOutlier(stats$subsets_Mito_percent, nmads=3, type="higher")
table(high_mito)
    # high_mito
    # FALSE  TRUE
    # 17245  2557

# Notes (from Erik):
# 1) Compare (location=="MT") vs grep("^MT-")
# 2) Drop based on qc metrics
# 3) Drop any genes with 0 counts in any nuclei 

# Getting cluster annotations from Erik's sce object
load(here("processed-data", "08_snRNA-seq_Erik", "s3e_hb.rda")) 

annoData <- Data.frame(row.names = colnames(s3e.hb), "SampleID" = s3e.hb$Sample,
                         colData(s3e.hb)$S




# 