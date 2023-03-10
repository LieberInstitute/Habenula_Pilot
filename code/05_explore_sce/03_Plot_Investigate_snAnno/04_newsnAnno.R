## March 8, 2023 - Bukola Ajanaku
# Updating snAnno and replotting to verify identities
# qrsh -l mem_free=20G,h_vmem=20G

library("SingleCellExperiment")
library("here")
library("tidyverse")
library("DeconvoBuddies")
library("sessioninfo")

# loading original sce object with just snAnno and other minor annotations
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
          "sce_post_09_clustered_qc.Rdata"))

# creating snAnno2 which will contain the combined MHb groups
sce$snAnno2 <- sce$snAnno

  ## LHb.6 is actually Endothelial. Total LHb is now 7 from 8.
sce$snAnno2[sce$snAnno2 == "LHb.6"] <- "Endo"
table(sce$snAnno2)
    # Astrocyte       Endo Excit.Thal         Hb Inhib.Thal      LHb.1      LHb.2 
    # 538         38       1800         51       7612        201        266 
    # LHb.3      LHb.4      LHb.5      LHb.7      LHb.8      MHb.1      MHb.2 
    # 134        477         83         39       1014        152        145 
    # MHb.3      MHb.4  Microglia      Oligo        OPC 
    # 395         18        145       2178       1796 

# doing the same for the  splitSNType. Has same info as snAnno but split by cluster (so annotation + cluster number)
sce$splitSNType2 <- sce$splitSNType
sce$splitSNType2[sce$splitSNType2 == "LHb.6_32"] <- "Endo"
table(sce$splitSNType2)
    # Astrocyte_11  Astrocyte_14  Astrocyte_34          Endo Excit.Thal_15 
    # 188           325            25            38           164 
    # Excit.Thal_18 Excit.Thal_19 Excit.Thal_20 Excit.Thal_21 Excit.Thal_22 
    # 373            51            86           217           276 
    # Excit.Thal_25 Excit.Thal_28 Excit.Thal_37  Excit.Thal_5  Excit.Thal_8 
    # 62            65            17           312           177 
    # Hb_30 Inhib.Thal_10 Inhib.Thal_17 Inhib.Thal_29 Inhib.Thal_35 
    # 51          2231          5241            85            17 
    # Inhib.Thal_6      LHb.1_13      LHb.2_24      LHb.3_26       LHb.4_3 
    # 38           201           266           134           477 
    # LHb.5_31      LHb.7_33       LHb.8_4      MHb.1_12      MHb.2_16 
    # 83            39          1014           152           145 
    # MHb.3_9      MHb.4_36   Microglia_1      Oligo_23      Oligo_27 
    # 395            18           145           280            65 
    # Oligo_7         OPC_2 
    # 1833          1796 

# # sourcing code from DLPFC Project (by Louise Huuki) 
   source(here("code", "04_snRNA-seq", "sourcing", "custom_plotExpression.R"))
# # actually runs and plots markers 
   source(here("code", "04_snRNA-seq", "sourcing", "my_plotMarkers.R"))

## new list of gene markers :)
# Terminal 1
new_markers.custom <- list(
  'Neuron' = c('SYT1', 'SNAP25', "SYT4", "SYP"), 
  "Non-Neuronal Subtype" = c("TNF", "KIR4.1", "KCNJ10"),
  'excitatory_neuron' = c('SLC17A6', 'SLC17A7'), 
  'inhibitory_neuron' = c('GAD1', 'GAD2'), 
  "Habenula_neurons" = c("CCBP2","CD63", "HTR5B", "KCNNH8", "KCTD8", "LRRC55", "MAPK4", "NEUROD1", 
                         "PIXNC1", "SCUBE1", "SSTR4", "TACR1", "SSTR2", "IRX2"),
  "LHB_neuron_specific" = c("PCDH10", "GABRA1", "SYN2", "GAP43", "HTR2C", "ADCYAP1R1", 
                            "CHRM3", "VGF", "GPR151", "SST", "SLC17AR"),
  "MHb_neuron_specific" = c("TAC3", "SPON1", "SEMA3d", "CALB1",'HHTR5b', "CNR1", "GPR4", "CHRNB3", "B4",
                            "SLC17A7"),
  'mediodorsal thalamus'= c('EPHA4','PDYN', 'LYPD6B', 'LYPD6', 'S1PR1', 'GBX2', 
                            'RAMP3', 'COX6A2', 'SLITRK6', 'DGAT2', "ADARB2"),
  "Endo/Mural" = c("CLDN5", "CARMN", "ITIH5", "NOTCH3", "ATP10A", "MECOM", "EBF1", 
                   "AC092957.1", "ITGA1", "VWF"),
  'oligodendrocyte' = c('MOBP', 'MBP', "CX3CR1"),
  'oligodendrocyte_precursor' = c('PDGFRA', 'VCAN'), 
  'microglia' = c('C3', 'CSF1R'), 
  'astrocyte' = c('GFAP', 'AQP4')
)
#  "Choroid Plexus" = c("klotho", "CLIC6", "OATP14", "EZRIN"),

extra_markers.custom <- list(
  "Macrophages" = c("CD14", "CD16", "CD64", "CD68", "CD71", "CCR5"),
  "Fibroblasts" = c("PDGFRA", "COL3A1", "Col1a1", "Col1a2", "Col5a1", "Loxl1", "Lum", "Fbln1", "Fbln2"),
  "Ependymal" = c("SLC6A11", "HDC", "Foxj1", "Pifo", "Dynlrb2"),
  "Pericytes" = c("ABCC9", "PDGFRB", "Cspg4"),
  "Polydendro" =  c("GPR17", "OLIG1", 'GAP43', 'PDGFRA')
)

eLife_markers.custom <- list(
  "MHb_Ventral_2thirds" = toupper(c("Fgf1", "Satb1", "Igfbp7", "Lmo3", "Slc18a3", "Tcf4", "Esam", "Chrna3", "Chrnb3")),
  "MHb_Ventrolateral" = c("Igfbp7", "Lmo3", "Slc18a3", "Tcf4", "Esam", "Syt15"),
  "MHb_Lateral" = c("Syt15", "Spon1", "Sema3d", "Calb1", "Rprm"),
  "MHb_Dorsal" = c("Tac2", "Calb1", "Rprm", "Col15a1", "Rasd2", "Adcyap1", "Wif1", "Cck", "Avail", "Fabp5"),
  "MHb_Superior" = c("Tac2", "Fgf1", "Satb1", "Col16a1", "Rasd2", "Adcyap1", "Wif1", "Cck", "Fxyd7", "Avil", 
                     "Asic4", "Pygm", "Fabp5", "Tac1"),
  "LHb_DEGs" = c("Gap43", "Rbfox1", "Parm1", "Chrnb3", "Bcl11b", "Th"),
  "LHb_HighlyExpressed_PerCluster" = c("Chrm3", "Vgf", "Gpr151", "Sst"),
  "LHb_Oval_Medial" = c("Gap43", "Gpd2", "Rbfox1", "Chrm3", "Vgf"),
  "LHb_Margina" = c("Gap43", 'Chrm3', "Paqr8", "Plch1", "Vgf", "Parm1"),
  "LHb_Lateral" = c("Pbx3", "Peg10", 'Parm1', "Gpr151", "Sst", "Cartpt"),
  "LHB_HbX" = c("Gpd2", "Gpr151", "Sst", "Cartpt", "Chrnb3", 'Bcl11b', "Th")
)

# Making all marker lists into all capital for plotexpression function below
new_markers.custom <- lapply(new_markers.custom, FUN = toupper)
extra_markers.custom <- lapply(extra_markers.custom, FUN = toupper)
eLife_markers.custom <- lapply(eLife_markers.custom, FUN = toupper)

# changing rownames for gene annotation purposes 
rownames(sce) <- rowData(sce)$Symbol

####### PLOTTING ###############################################################
# creating plot dir
plot_dir <- here("plots", "05_explore_sce", "03_Plot_Investigate_snAnno", 
                 "04_newsnAnno")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

#### For snAnno2 (updated snAnno) ####
# 1) New_markers.custom 
  snAnnoCustom_new <- here(plot_dir, "snAnno_new_custom_markers_violin_plots.pdf")
  my_plotMarkers(sce, marker_list = new_markers.custom, assay = "logcounts", 
                 cat = "snAnno2", fill_colors = NULL, pdf_fn = snAnnoCustom_new)
   # test <-  plot_marker_express_List(sce, 
   #                           pdf_fn =  snAnnoCustom_new, 
   #                           gene_list = new_markers.custom,
   #                           cat = "snAnno2") 
 
# 2) Extra Markers doesn't work for snAnno combination   
      
# 3) eLife marker categories  
  eLife_mark <- here(plot_dir, "snAnno_eLife_categories_violin_plots.pdf")
    my_plotMarkers(sce, marker_list = eLife_markers.custom, assay = "logcounts", 
                   cat = "snAnno2", fill_colors = NULL, pdf_fn = eLife_mark)

#### For splitSNType2 (updated splitSNType) ####
# 1) new marker list on all 37 clusters to double check higher order combos
snAnnoCustom_newer1 <- here(plot_dir, "SPLIT_snAnno_new_custom_markers_violin_plots.pdf")
my_plotMarkers(sce, marker_list = new_markers.custom, assay = "logcounts", 
                   cat = "splitSNType2", fill_colors = NULL, pdf_fn = snAnnoCustom_newer1)

# 2) extra marker list to make sure we aren't missing any smalller classes for higher 
# order combination
extraMarkers1 <- here(plot_dir, "SPLIT_snAnno_extra_marker_categories_violin_plots.pdf")
my_plotMarkers(sce, marker_list = extra_markers.custom, assay = "logcounts", 
               cat = "splitSNType2", fill_colors = NULL, pdf_fn = extraMarkers1)

# 3) eLife marker categorieS just for fun  
eLife_mark_newer1 <- here(plot_dir, "SPLIT_snAnno_eLife_categories_violin_plots.pdf")
my_plotMarkers(sce, marker_list = eLife_markers.custom, assay = "logcounts", 
               cat = "splitSNType2", fill_colors = NULL, pdf_fn = eLife_mark_newer1)
# Saving 
save(sce, file = here("processed-data", "04_snRNA-seq", "sce_objects", "sce_with_snAnno2.RDATA"))


# # updating for combined MHb.3 and MHb.2 for snAnno3 
# sce$snAnno3 <- sce$snAnno2
# # combining MHb.3 with MHb.2
# sce$snAnno3[sce$snAnno3 == "MHb.3"] <- "MHb.2"
# 
# sce$splitSNType3 <- sce$splitSNType2
# sce$splitSNType3[sce$splitSNType3 == "MHb.3_9"] <- "MHb.2_16"
# save(sce, file = here("processed-data", "04_snRNA-seq", "sce_objects", "sce_with_snAnno2_and_snAnno3.RDATA"))


