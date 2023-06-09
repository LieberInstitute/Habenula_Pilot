## June 9, 2023 - Bukola Ajanaku
# Calculating t-stats for trans-special (mouse vs humann) habenula data.
# qrsh -l mem_free=50G,h_vmem=50G

library("SingleCellExperiment")
library("here")
library("sessioninfo")
library("spatialLIBD")
# library("Seurat")

read.csv(file = here("processed-data",  
                 "09_trans_special_analysis",
                 "habData.csv"))


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

sum(startsWith(colnames(sce.mm.sub), "hab_161105"))

