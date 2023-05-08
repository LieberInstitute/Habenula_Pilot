## May 8, 2023 - Bukola Ajanaku
# Fixing pre-drop clusters by encorporating OPC_noidy as an annotated
# class.
# qrsh -l mem_free=50G,h_vmem=50G

library("SingleCellExperiment")
library("here")
library("tidyverse")
library("scater")
library("sessioninfo")

# loading sce object with non dropped final annotations
load(file = here("processed-data", "99_paper_figs", "sce_objects", 
                 "sce_final_preHbdrop.RDATA"))
sce <- sce_final_preHbdrop
rm(sce_final_preHbdrop)

dim(sce)
  # [1] 33848 17082

table(sce$final_Annotations)
  # Astrocyte    Endo Excit.Neuron   Excit.Thal   Inhib.Thal        LHb.1 
  # 538           38           51         1800         7612          201 
  # LHb.2        LHb.3        LHb.4        LHb.5        LHb.6        LHb.7 
  # 266          134          477           83           39         1014 
  # MHb.1        MHb.2        MHb.3    Microglia        Oligo          OPC 
  # 152          540           18          145         2178         1796 

# changing OPC_noisy class into a cluster of it's own.
# adding rownames of colData as a row for easier subsetting
OPC_noisy_Samps = c("Br5555", "Br1204", "Br1092")

# grabbing barcodes for  OPC
onlyOPC <- sce[, which(sce$final_Annotations == "OPC")]

# grabbing noisy OPC
onlyOPC <- onlyOPC[, which(onlyOPC$Sample %in% OPC_noisy_Samps)]

# grabbing barcodes            
RowNos <- rownames(colData(onlyOPC))

# creating OPC_noisy subset
sce[, rownames(colData(sce)) %in% RowNos]$final_Annotations <- "OPC_noisy"

# check
  # table(sce$Sample, sce$final_Annotations)
  # Astrocyte Endo Excit.Neuron Excit.Thal Inhib.Thal LHb.1 LHb.2 LHb.3
  # Br1092         4    0            4        959       1716    19    25    17
  # Br1204         7    1           46        194        117    35    69    44
  # Br1469       237    4            1          5          1     1     0     0
  # Br1735       151   18            0         51       1825    10    25    16
  # Br5555         3    0            0        515       1181   128   147    57
  # Br5558         0    0            0         32        746     0     0     0
  # Br5639       136   15            0         44       2026     8     0     0
  # 
  # LHb.4 LHb.5 LHb.6 LHb.7 MHb.1 MHb.2 MHb.3 Microglia Oligo  OPC
  # Br1092    27    15     0    34    25    27     2         2    13    0
  # Br1204   127     9    13   275    96    61     2         0     9    0
  # Br1469     0     0     0    14     0    16     0       107  1435  330
  # Br1735    47    14     0   133     2     8     3        15   351  432
  # Br5555   275    45    26   555    29   425    11         1     7    0
  # Br5558     0     0     0     0     0     0     0         0     0    0
  # Br5639     1     0     0     3     0     3     0        20   363  440
  # 
  # OPC_noisy
  # Br1092       323
  # Br1204       134
  # Br1469         0
  # Br1735         0
  # Br5555       137
  # Br5558         0
  # Br5639         0


# save
save(sce, file = here("processed-data", "99_paper_figs", "sce_objects", 
                      "sce_final_preHbdrop.RDATA"))

# 