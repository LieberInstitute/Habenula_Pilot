## March 8, 2023 - Bukola Ajanaku
# Updating snAnno and replotting to verify identities
# qrsh -l mem_free=20G,h_vmem=20G

library("SingleCellExperiment")
library("here")
library("tidyverse")

# loading original sce object with just snAnno and other minor annotations
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
          "sce_post_09_clustered_qc.Rdata"))

# creating snAnno2 which will contain the combined MHb groups
sce$snAnno2 <- sce$snAnno

  ## LHb.6 is actually Endothelial. Total LHb is now 7 from 8.
sce$snAnno2[sce$snAnno2 == "LHB.6"] <- "Endo"

# sourcing code from DLPFC Project (by Louise Huuki) 
  source(here("code", "04_snRNA-seq", "sourcing", "custom_plotExpression.R"))
# actually runs and plots markers 
  source(here("code", "04_snRNA-seq", "sourcing", "my_plotMarkers.R"))

## new list of gene markers :)
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
  "Choroid Plexus" = c("klotho", "CLIC6", "OATP14", "EZRIN"),
  'oligodendrocyte' = c('MOBP', 'MBP', "CX3CR1"),
  'oligodendrocyte_precursor' = c('PDGFRA', 'VCAN'), 
  'microglia' = c('C3', 'CSF1R'), 
  'astrocyte' = c('GFAP', 'AQP4')
)

extra_markers.custom <- list(
  "Macrophages" = c("MFC1", "MFC1"),
  "Fibroblasts" = c("PDGFRA", "COL3A1"),
  "Ependymal" = c("SLC6A11", "HDC"),
  "Pericytes" = c("ABCC9", "PDGFRB"),
  "Polydendro" =  c("GPR17", "OLIGO1", 'GAP43', 'PDGFRA')
)


eLife_markers.custom <- list(
  "MHb_Ventral_2thirds" = c("Fgf1", "Satb1", "Igfbp7", "Lmo3", "Slc18a3", "Tcf4", "Esam", "Chrna3", "Chrnb3"),
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

# changing rownames for gene annotation purposes 
rownames(sce) <- rowData(sce)$Symbol


####### PLOTTING ###############################################################
# creating plot dir
plot_dir <- here("plots", "05_explore_sce", "03_Plot_Investigate_snAnno", 
                 "04_newsnAnno")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

# 1) overall new_markers.custom creating name of pdf
  snAnnoCustom_new <- here(plot_dir, "snAnno_new_custom_markers_violin_plots.pdf")
  # plotting by snAnno in sce because snAnno is my original combined cluster resolution for sn annotations
    my_plotMarkers(sce, marker_list = new_markers.custom, assay = "logcounts", 
                   cat = "snAnno2", fill_colors = NULL, pdf_fn = snAnnoCustom_new)
    
# 2) throwing in some extra marker categories 
  extraMark <- here(plot_dir, "snAnno_extra_marker_categories_violin_plots.pdf")
  # plotting by snAnno in sce because snAnno is my original combined cluster resolution for sn annotations
      my_plotMarkers(sce, marker_list = new_markers.custom, assay = "logcounts", 
                     cat = "snAnno2", fill_colors = NULL, pdf_fn = extraMark)
      
# 3) eLife marker categories 
  eLife_mark <- here(plot_dir, "snAnno_eLife_categories_violin_plots.pdf")
# plotting by snAnno in sce because snAnno is my original combined cluster resolution for sn annotations
    my_plotMarkers(sce, marker_list = new_markers.custom, assay = "logcounts", 
                   cat = "snAnno2", fill_colors = NULL, pdf_fn = eLife_mark)


# tell Louise that sometimes code can get hooked when missing markers. Fix is to cancel function and then dev.off()


