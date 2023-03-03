## March 3, 2023 - Bukola Ajanaku
# Running mean ratio expression analysis to help confirm annotated identity of 
# kept clusters (all 37).
# qrsh -l mem_free=20G,h_vmem=20G

library("SingleCellExperiment")
library("here")

# loading sce with annotation names 
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
          "sce_post_09_clustered_qc.Rdata"))

dim(sce)
    # [1] 33848 17082

table(sce$Sample)
    # Br1092 Br1204 Br1469 Br1735 Br5555 Br5558 Br5639 
    # 3212   1239   2151   3101   3542    778   3059 

# We are using WalkTrap 10
table(sce$wT_10_Erik)
    # 10wTrap_1 10wTrap_10 10wTrap_11 10wTrap_12 10wTrap_13 10wTrap_14 10wTrap_15 
    # 145       2231        188        152        201        325        164 
    # 10wTrap_16 10wTrap_17 10wTrap_18 10wTrap_19  10wTrap_2 10wTrap_20 10wTrap_21 
    # 145       5241        373         51       1796         86        217 
    # 10wTrap_22 10wTrap_23 10wTrap_24 10wTrap_25 10wTrap_26 10wTrap_27 10wTrap_28 
    # 276        280        266         62        134         65         65 
    # 10wTrap_29  10wTrap_3 10wTrap_30 10wTrap_31 10wTrap_32 10wTrap_33 10wTrap_34 
    # 85        477         51         83         38         39         25 
    # 10wTrap_35 10wTrap_36 10wTrap_37  10wTrap_4  10wTrap_5  10wTrap_6  10wTrap_7 
    # 17         18         17       1014        312         38       1833 
    # 10wTrap_8  10wTrap_9 
    # 177        395 

# Each cluster was annotated as so:
table(sce$splitSNType)
    # Astrocyte_11  Astrocyte_14  Astrocyte_34 Excit.Thal_15 Excit.Thal_18 
    # 188           325            25           164           373 
    # Excit.Thal_19 Excit.Thal_20 Excit.Thal_21 Excit.Thal_22 Excit.Thal_25 
    # 51            86           217           276            62 
    # Excit.Thal_28 Excit.Thal_37  Excit.Thal_5  Excit.Thal_8         Hb_30 
    # 65            17           312           177            51 
    # Inhib.Thal_10 Inhib.Thal_17 Inhib.Thal_29 Inhib.Thal_35  Inhib.Thal_6 
    # 2231          5241            85            17            38 
    # LHb.1_13      LHb.2_24      LHb.3_26       LHb.4_3      LHb.5_31 
    # 201           266           134           477            83 
    # LHb.6_32      LHb.7_33       LHb.8_4      MHb.1_12      MHb.2_16 
    # 38            39          1014           152           145 
    # MHb.3_9      MHb.4_36   Microglia_1      Oligo_23      Oligo_27 
    # 395            18           145           280            65 
    # Oligo_7         OPC_2 
    # 1833          1796 

# In summary:
table(sce$snAnno)

  # Astrocyte Excit.Thal         Hb Inhib.Thal      LHb.1      LHb.2      LHb.3 
  # 538       1800         51       7612        201        266        134 
  # LHb.4      LHb.5      LHb.6      LHb.7      LHb.8      MHb.1      MHb.2 
  # 477         83         38         39       1014        152        145 
  # MHb.3      MHb.4  Microglia      Oligo        OPC 
  # 395         18        145       2178       1796 