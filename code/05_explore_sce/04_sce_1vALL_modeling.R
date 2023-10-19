
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

table(sce$final_Annotations)
# Astrocyte       Endo Excit.Thal Inhib.Thal      LHb.1      LHb.2      LHb.3      LHb.4      LHb.5      LHb.6 
#       538         38       1800       7612        201        266        134        477         83         39 
# LHb.7      MHb.1      MHb.2      MHb.3  Microglia      Oligo        OPC 
# 1014        152        540         18        145       2178       1202 

table(sce$broad_Annotations)
# Astrocyte Endo Excit.Thal Inhib.Thal        LHb        MHb  Microglia      Oligo        OPC 
# 538         38       1800       7612       2214        710        145       2178       1202 

## Use ENSEMBL IDs as rownames
rownames(sce) <- rowData(sce)$ID

## run for final_Annotations
message(Sys.time(), "- Run registration_wrapper final_Annotations")
sce_modeling_final_Annotations <- registration_wrapper(
  sce,
  var_registration = "final_Annotations",
  var_sample_id = "Sample",
  covars = NULL,
  gene_ensembl = "ID",
  gene_name = "Symbol",
  suffix = "",
  min_ncells = 10,
  pseudobulk_rds_file = here(data.dir, "sce_pseudo_final_Annotations.rds")
)

save(sce_modeling_final_Annotations, file = here(data.dir, "sce_modeling_final_Annotations.Rdata"))


## run for broad_Annotations
message(Sys.time(), "- Run registration_wrapper broad_Annotations")
sce_modeling_broad_Annotations <- registration_wrapper(
  sce,
  var_registration = "broad_Annotations",
  var_sample_id = "Sample",
  covars = NULL,
  gene_ensembl = "ID",
  gene_name = "Symbol",
  suffix = "",
  min_ncells = 10,
  pseudobulk_rds_file = here(data.dir, "sce_pseudo_broad_Annotations.rds")
)

save(sce_modeling_broad_Annotations, file = here(data.dir, "sce_modeling_broad_Annotations.Rdata"))

#### Broad 2 - combine Habenula ####

sce$broad_Annotations2 <- gsub("^.*Hb", "Hb",sce$broad_Annotations)
table(sce$broad_Annotations2)
# Astrocyte       Endo Excit.Thal         Hb Inhib.Thal  Microglia      Oligo        OPC 
#       538         38       1800       2924       7612        145       2178       1202 

message(Sys.time(), "- Run registration_wrapper broad_Annotations2")
sce_modeling_broad_Annotations2 <- registration_wrapper(
  sce,
  var_registration = "broad_Annotations2",
  var_sample_id = "Sample",
  covars = NULL,
  gene_ensembl = "ID",
  gene_name = "Symbol",
  suffix = "",
  min_ncells = 10,
  pseudobulk_rds_file = here(data.dir, "sce_pseudo_broad_Annotations2.rds")
)

save(sce_modeling_broad_Annotations2, file = here(data.dir, "sce_modeling_broad_Annotations2.Rdata"))


# sgejobs::job_single('04_sce_1vALL_modeling', create_shell = TRUE, memory = '30G', command = "Rscript 04_sce_1vALL_modeling.R")

## Reproducibility information
print("Reproducibility information:")
options(width = 120)
session_info()
