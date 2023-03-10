## 2/21/23 - Bukola Ajanaku
# Pseudobulking to on snAnno level to investigate new annotations 
# qrsh -l mem_free=20G,h_vmem=20G
# Terminal 7

library("SingleCellExperiment")
library("here")
library("spatialLIBD")

load(here("processed-data", "04_snRNA-seq", "sce_objects", 
          "sce_post_09_clustered_qc.Rdata"))

# updating LHb.6 as Endo 
# creating snAnno2 which will contain the combined MHb groups
sce$snAnno2 <- sce$snAnno

## LHb.6 is actually Endothelial. Total LHb is now 7 from 8.
sce$snAnno2[sce$snAnno2 == "LHb.6"] <- "Endo"

## Combining MHb.3 with MHb.2
sce$snAnno3 <- sce$snAnno2

# combining MHb.3 with MHb.2
sce$snAnno3[sce$snAnno3 == "MHb.3"] <- "MHb.2"


# Pseudobulking to create compressed sce object
message("Start - snAnno2")
Sys.time()

  sce_new_pb_snAnno2 <- registration_pseudobulk(sce, "snAnno2", "Sample")

message("End - Pseudobulk Complete")
Sys.time()

message("Start - snAnno3")
Sys.time()

  sce_new_pb_snAnno3 <- registration_pseudobulk(sce, "snAnno3", "Sample")

message("End - Pseudobulk Complete")
Sys.time()

save(sce_new_pb_snAnno2, sce_new_pb_snAnno3, file = here("processed-data", "04_snRNA-seq", "sce_objects", 
                                            "sce_new_pseudobulk_with_snAnno_version2and3.Rdata"))

# sgejobs::job_single('07c_Pseudobulking', create_shell = TRUE, queue= 'bluejay', memory = '25G', command = "Rscript 07c_Pseudobulking.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


