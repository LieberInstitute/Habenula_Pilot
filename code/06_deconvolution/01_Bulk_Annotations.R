## March 17, 2023 - Bukola Ajanaku
# Here I will be creating different collapse combinations for final_Annotations (the official
# clustering method we are using on my cleaned sce object) for my bulk Annotations!
# Based on my work in 05_Updated_Annotations_meanExpression.R
# qrsh -l mem_free=20G,h_vmem=20G

library("SingleCellExperiment")
library("here")

# loading regular sce object
load(here("processed-data", "99_paper_figs", "sce_objects", 
          "official_final_sce.RDATA"), verbose = TRUE)

##### Cleaning up sce object ###################################################
## check 
table(sce$final_Annotations)
    # Astrocyte   Endo Excit.Thal Inhib.Thal      LHb.1      LHb.2      LHb.3 
    # 538         38       1800       7612        201        266        134 
    # LHb.4       LHb.5      LHb.6      LHb.7      MHb.1      MHb.2      MHb.3 
    # 477         83         39       1014        152        540         18 
    # Microglia  Oligo        OPC 
    # 145        2178       1796

table(sce$OPC_clean)
# No   Yes 
# 594 16437

# dropping noisy OPC
sce <- sce[, which(sce$OPC_clean == "Yes")]

table(sce$OPC_clean)
# Yes 
# 16437 

##### Bulk Annotations #########################################################
# 1: bulkTypeSepHb 
# Combining into 5 glia, two thalamus, 1 broad LHb, and 1 broad MHb.
sce$bulkTypeSepHb <- sce$final_Annotations
  # making separated Hb (2)
sce$bulkTypeSepHb[sce$bulkTypeSepHb %in% grep("^LHb\\.", unique(sce$bulkTypeSepHb), value = TRUE)] <-  'LHb'
sce$bulkTypeSepHb[sce$bulkTypeSepHb %in% grep("^MHb\\.", unique(sce$bulkTypeSepHb), value = TRUE)] <-  'MHb'

  # check!
  table(sce$bulkTypeSepHb)
    # Astrocyte       Endo Excit.Thal Inhib.Thal        LHb        MHb  Microglia 
    # 538         38       1800       7612       2214        710        145 
    # Oligo        OPC 
    # 2178       1796


# 2: bulkTypeBroadHb
# Combining into 5 glia, two thalamus, 1 broad Hb.
sce$bulkTypeBroadHb <- sce$final_Annotations
  # making one big Hb
sce$bulkTypeBroadHb[sce$bulkTypeBroadHb %in% grep("^LHb\\.", unique(sce$bulkTypeBroadHb), value = TRUE)] <-  'Hb'
sce$bulkTypeBroadHb[sce$bulkTypeBroadHb %in% grep("^MHb\\.", unique(sce$bulkTypeBroadHb), value = TRUE)] <-  'Hb'

  # check!
  table(sce$bulkTypeBroadHb)
    # Astrocyte   Endo.  Excit.Thal    Hb   Inhib.Thal  Microglia      Oligo 
    # 538         38       1800       2924       7612        145       2178 
    # OPC 
    # 1796

# 3: bulkTypeAllCollapse
# Combining into 5 glia, one thalamus, 1 broad Hb.
sce$bulkTypeAllCollapse <- sce$final_Annotations
  # making one big Hb
sce$bulkTypeAllCollapse[sce$bulkTypeAllCollapse %in% grep("^LHb\\.", unique(sce$bulkTypeAllCollapse), value = TRUE)] <-  'Hb'
sce$bulkTypeAllCollapse[sce$bulkTypeAllCollapse %in% grep("^MHb\\.", unique(sce$bulkTypeAllCollapse), value = TRUE)] <-  'Hb'
  # combining thals
sce$bulkTypeAllCollapse[sce$bulkTypeAllCollapse %in% grep("^Excit\\.", unique(sce$bulkTypeAllCollapse), value = TRUE)] <-  'Thal'
sce$bulkTypeAllCollapse[sce$bulkTypeAllCollapse %in% grep("^Inhib\\.", unique(sce$bulkTypeAllCollapse), value = TRUE)] <-  'Thal'

  # check! 
  table(sce$bulkTypeAllCollapse)
    # Astrocyte   Endo    Hb        Microglia  Oligo       OPC      Thal 
    # 538         38      2924       145        2178      1202      9412 

new_dir <- here("processed-data", "06_deconvolution", "sce_objects")
if(!dir.exists(new_dir)){
  dir.create(new_dir) # recursive = TRUE
}

# saving sce object
save(sce, file = here(new_dir, "sce_first_bulkTypes.RDATA"))

