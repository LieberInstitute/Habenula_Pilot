## 2/21/23 - Bukola Ajanaku
# Pseudobulking to help with gene annotations
# qrsh -l mem_free=20G,h_vmem=20G

library("SingleCellExperiment")
library("here")
library("spatialLIBD")

load(here("processed-data", "04_snRNA-seq", "sce_objects", 
          "sce_post_clustering_with_logcounts.Rdata"))

# Pseudobulking to create compressed sce object

# WalkTrap 50
message("Start - Pseudobulk 10")
  Sys.time()
sce_psuedo_wT10 <- registration_pseudobulk(sce, "wT_10_Erik", "Sample")

message("Start - Pseudobulk 20")
  Sys.time()
sce_psuedo_wT20 <- registration_pseudobulk(sce, "wT_20_Erik", "Sample")
  
message("Start - Pseudobulk 50")
  Sys.time()
sce_psuedo_wT50 <- registration_pseudobulk(sce, "wT_50_Erik", "Sample")

message("End - Pseudobulk Complete")
  Sys.time()
  
save(sce_psuedo_wT10, sce_psuedo_wT20, sce_psuedo_wT50, file = here("processed-data", "04_snRNA-seq", "sce_objects", 
                                                                    "sce_pseduobulked_wT.Rdata"))
     
# sgejobs::job_single('07c_Pseudobulking', create_shell = TRUE, queue= 'bluejay', memory = '25G', command = "Rscript 07c_Pseudobulking.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


