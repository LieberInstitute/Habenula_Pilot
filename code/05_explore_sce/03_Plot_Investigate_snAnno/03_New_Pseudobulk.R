## 2/21/23 - Bukola Ajanaku
# Pseudobulking to on snAnno level to investigate new annotations 
# qrsh -l mem_free=20G,h_vmem=20G

library("SingleCellExperiment")
library("here")
library("spatialLIBD")

load(here("processed-data", "04_snRNA-seq", "sce_objects", 
          "sce_post_09_clustered_qc.Rdata"))

# Pseudobulking to create compressed sce object
message("Start - snAnno")
Sys.time()

  sce_new_pb_snAnno <- registration_pseudobulk(sce, "snAnno", "Sample")

message("End - Pseudobulk Complete")
Sys.time()


save(sce_new_pb_snAnno, file = here("processed-data", "04_snRNA-seq", "sce_objects", 
                                            "sce_new_pseudobulk_with_snAnno.Rdata"))

# sgejobs::job_single('07c_Pseudobulking', create_shell = TRUE, queue= 'bluejay', memory = '25G', command = "Rscript 07c_Pseudobulking.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


