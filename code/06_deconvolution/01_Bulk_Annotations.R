## March 17, 2023 - Bukola Ajanaku
# Here I will be creating different collapse combinations for snAnno3 (the official
# clustering method we are using on my cleaned sce object) for my bulk Annotations!
# Based on my work in 05_Updated_Annotations_meanExpression.R
# qrsh -l mem_free=20G,h_vmem=20G

library("SingleCellExperiment")
library("here")

# loading regular sce object
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
          "sce_with_snAnno2_and_snAnno3.RDATA"))

##### Cleaning up sce object ###################################################
dropNames <- c("louvain10", "louvain20", "louvain50", "wT_10_Erik", "wT_20_Erik",
               "wT_50_Erik", "groupErik", "wt10_ANNO", "wt50_ANNO", "cellType_wT10",
               "cellType_wT20", "cellType_wT50", "splitProbClusts", "splitSNType", 
               "snAnno", "snAnno2", "splitSNType2", "splitSNType3")

for(i in dropNames){
  colData(sce)[, dropNames] <- NULL
}

table(sce$snAnno3)
  # Astrocyte       Endo Excit.Thal         Hb Inhib.Thal      LHb.1      LHb.2 
  # 538         38       1800         51       7612        201        266 
  # LHb.3      LHb.4      LHb.5      LHb.7      LHb.8      MHb.1      MHb.2 
  # 134        477         83         39       1014        152        540 
  # MHb.4  Microglia      Oligo        OPC 
  # 18        145       2178       1796

# Shifting columns upward where there was a vacancy. Not sure why this sce object 
# does not already incude this shift made earier.
sce$snAnno3[sce$snAnno3 == "LHb.7"] <- "LHb.6"
sce$snAnno3[sce$snAnno3 == "LHb.8"] <- "LHb.7"
sce$snAnno3[sce$snAnno3 == "MHb.4"] <- "MHb.3"

table(sce$snAnno3)
  # Astrocyte       Endo Excit.Thal         Hb Inhib.Thal      LHb.1      LHb.2 
  # 538         38       1800         51       7612        201        266 
  # LHb.3      LHb.4      LHb.5      LHb.6      LHb.7      MHb.1      MHb.2 
  # 134        477         83         39       1014        152        540 
  # MHb.3  Microglia      Oligo        OPC 
  # 18        145       2178       1796 

# saving sce object
save(sce, file = here("processed-data", "04_snRNA-seq", "sce_objects", "final_sce_official_annotation.RDATA"))

##### Bulk Annotations #########################################################
# bulkAnno1 
# Combining into 








#



