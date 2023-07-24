
library("SingleCellExperiment")
library("spatialLIBD")
library("here")
library("sessioninfo")

data.dir <- here("processed-data", "05_explore_sce","04_sce_1vALL_modeling")
if(!dir.exists(data.dir)) dir.create(data.dir)

load(here("processed-data", "sce_objects", "sce_Habenula_Pilot.Rdata"), verbose = TRUE)
#sce

dim(sce)
# [1] 33848 16437

colnames(colData(sce))

table(sce$final_Annotations)
# Astrocyte       Endo Excit.Thal Inhib.Thal      LHb.1      LHb.2      LHb.3      LHb.4      LHb.5      LHb.6 
#       538         38       1800       7612        201        266        134        477         83         39 
# LHb.7      MHb.1      MHb.2      MHb.3  Microglia      Oligo        OPC 
# 1014        152        540         18        145       2178       1202 

table(sce$final_Annotations)

message(Sys.time(), "- Run registration_wrapper")

sce_modeling_results <- registration_wrapper(
  sce,
  var_registration = "final_Annotations",
  var_sample_id = "Sample",
  covars = NULL,
  gene_ensembl = NULL,
  gene_name = NULL,
  suffix = "",
  min_ncells = 10,
  pseudobulk_rds_file = here(data.dir, "sce_pseudo_final_Annotations.rds")
)

save(sce_modeling_results, file = here(data.dir, "sce_modeling_results.Rdata"))

# sgejobs::job_single('04_sce_1vALL_modeling', create_shell = TRUE, memory = '30G', command = "Rscript 04_sce_1vALL_modeling.R")

## Reproducibility information
print("Reproducibility information:")
options(width = 120)
session_info()
