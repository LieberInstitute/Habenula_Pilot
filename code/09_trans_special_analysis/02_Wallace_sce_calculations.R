## June 9, 2023 - Bukola Ajanaku
# Calculating t-stats for trans-special (mouse vs humann) habenula data.
# qrsh -l mem_free=50G,h_vmem=50G

library("SingleCellExperiment")
library("here")
library("sessioninfo")
library("spatialLIBD")
# library("Seurat")

# loading the homologous sce objects
load(file = here("processed-data",  
                 "09_trans_special_analysis",
                 "sce_homologs_mm_hsap.rda"), verbose = TRUE)
  # sce.mm.sub
  # sce.hsap.sub

table(sce.hsap.sub$final_Annotations)
# Astrocyte       Endo Excit.Thal Inhib.Thal      LHb.1      LHb.2      LHb.3 
# 538         38       1800       7612        201        266        134 
# LHb.4      LHb.5      LHb.6      LHb.7      MHb.1      MHb.2      MHb.3 
# 477         83         39       1014        152        540         18 
# Microglia      Oligo        OPC 
# 145       2178       1202 

table(sce.mm.sub$orig.ident)
# hab 
# 7506 <- only habenula data in this data set

## NOTE: 
# query = mouse
# reference = human

# rownames(sce.mm.sub) <- rowData(sce.mm.sub)$JAX.geneID
# 
# sce_modeling_results_hsap <- registration_wrapper(
#   sce = sce.hsap.sub,
#   var_registration = "final_Annotations",
#   var_sample_id = "Sample",
#   gene_ensembl = "ID",
#   gene_name = "Symbol"
# )
# 
# sce_pseudo_mm <- registration_pseudobulk(
#    sce = sce.mm.sub,
#    var_registration = "ident",
#    var_sample_id = "Sample",
#    min_ncells = NULL
# )
# 
# 
# reg.mod.mm <- registration_model(
#   sce.mm.sub,
#   covars = NULL,
#   var_registration = "Sample"
# )












# .


