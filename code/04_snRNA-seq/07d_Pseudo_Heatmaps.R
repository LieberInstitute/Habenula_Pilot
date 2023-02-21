## 2/21/23 - Bukola Ajanaku
# Heatmapping pseudobulked data against marker gene list.
# qrsh -l mem_free=20G,h_vmem=20G

library("SingleCellExperiment")
library("here")
library("sessioninfo")
library("ComplexHeatmap")
library("tidyverse")

# must be pseudo_bulked data: 
## sce_psuedo_wT10, sce_psuedo_wT20, sce_psuedo_wT50
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
          "sce_pseduobulked_wT.Rdata"), verbose = TRUE)

# list of marker genes 
markers.custom = list(
  'neuron' = c('SYT1', 'SNAP25'), # 'GRIN1','MAP2'),
  'excitatory_neuron' = c('SLC17A6', 'SLC17A7'), # 'SLC17A8'),
  'inhibitory_neuron' = c('GAD1', 'GAD2'), #'SLC32A1'),
  'mediodorsal thalamus'= c('EPHA4','PDYN', 'LYPD6B', 'LYPD6', 'S1PR1', 'GBX2', 
                            'RAMP3', 'COX6A2', 'SLITRK6', 'DGAT2', "ADARB2"), # LYPD6B and ADARB2*, EPHA4 may not be best
  'Pf/PVT' = c("INHBA", "NPTXR"), 
  'Hb neuron specific'= c('POU4F1','GPR151'), #'STMN2','NR4A2','VAV2','LPAR1'), #CALB2 and POU2F2 were thrown out because not specific enough]
  'MHB neuron specific' = c('TAC1','CHAT','CHRNB4', "TAC3"),#'SLC17A7'
  'LHB neuron specific' = c('HTR2C','MMRN1', "ANO3"),#'RFTN1'
  'oligodendrocyte' = c('MOBP', 'MBP'), # 'PLP1'),
  'oligodendrocyte_precursor' = c('PDGFRA', 'VCAN'), # 'CSPG4', 'GPR17'),
  'microglia' = c('C3', 'CSF1R'), #'C3'),
  'astrocyte' = c('GFAP', 'AQP4'),
  "Endo/CP" = c("TTR", "FOLR1", "FLT1", "CLDN5")
)
# "FOLR1" does not have high enough expression

# Making data frame
markTable <- as.data.frame(unlist(markers.custom)) |> 
  rownames_to_column("cellType") |>
  rename(gene = `unlist(markers.custom)`) |>
  mutate(cellType = gsub("\\d+", "", cellType)) |>
  filter(gene %in% rowData(sce_psuedo_wT10)$Symbol)
# "FOLR1" does not have high enough expression in sce10

# Replacing genes with symbols for heatmap
rownames(sce_psuedo_wT10) <- rowData(sce_psuedo_wT10)$Symbol

# extracting logcounts from sce and subset just genes in markTable
markerlogcounts_wT10 <- logcounts(sce_psuedo_wT10[markTable$gene, ])

# getting z scores
marker_z_score <- scale(t(markerlogcounts_wT10))
# corner(marker_z_score)

# heatmap
plot_dir <- here("plots", "04_snRNA-seq", "07d_Pseudo_Heatmaps")
  if(!dir.exists(plot_dir)){
    dir.create(plot_dir)
  }

column_ha <- HeatmapAnnotation(
  cell_type = markTable$cellType
)

row_ha <- rowAnnotation(
  cluster = sce_psuedo_wT10$wT_10_Erik
)

pdf(here(plot_dir, "markers_heatmap_layer_wT10.pdf"), width = 17, height = 11)
  Heatmap(marker_z_score,
          cluster_rows = TRUE,
          cluster_columns = FALSE,
          right_annotation = row_ha,
          top_annotation = column_ha,
          column_split = markTable$cellType)
dev.off()





# 



