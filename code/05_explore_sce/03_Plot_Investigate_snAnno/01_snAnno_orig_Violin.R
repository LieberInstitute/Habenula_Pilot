## March 8, 2023 - Bukola Ajanaku
# Taking the originally combined clusters from walktrap10 with the resolution 
# necessary for signle nucleus annotation in the lateral and medial habenula and 
# violin plotting against different sets of markers in order to clearly validate their 
# their identities.
# qrsh -l mem_free=20G,h_vmem=20G

library("SingleCellExperiment")
library("here")
library("sessioninfo")
library("jaffelab")
library("ggplot2")

####### LOADING ################################################################
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
          "sce_post_09_clustered_qc.Rdata"))

####### SOURCING ###############################################################
# sourcing code from DLPFC Project (by Louise Huuki) 
source(here("code", "04_snRNA-seq", "sourcing", "custom_plotExpression.R"))

# actually runs and plots markers 
source(here("code", "04_snRNA-seq", "sourcing", "my_plotMarkers.R"))

####### BODY ###################################################################
# My custom marker genes list (heavily influenced by Erik Nelson and Louise Huuki)
markers.custom = list(
  'neuron' = c('SYT1', 'SNAP25'), 
  'excitatory_neuron' = c('SLC17A6', 'SLC17A7'), 
  'inhibitory_neuron' = c('GAD1', 'GAD2'), 
  'mediodorsal thalamus'= c('EPHA4','PDYN', 'LYPD6B', 'LYPD6', 'S1PR1', 'GBX2', 
                            'RAMP3', 'COX6A2', 'SLITRK6', 'DGAT2', "ADARB2"), 
  'Pf/PVT' = c("INHBA", "NPTXR"), 
  'Hb neuron specific'= c('POU4F1','GPR151', 'VAV2', "NR4A2", "NEUROD1"), 
  'MHB neuron specific' = c('TAC1','CHAT','CHRNB4', "TAC3"),
  'LHB neuron specific' = c('HTR2C','MMRN1', "ANO3", "ESRRG"),
  'oligodendrocyte' = c('MOBP', 'MBP'),
  'oligodendrocyte_precursor' = c('PDGFRA', 'VCAN'), 
  'microglia' = c('C3', 'CSF1R'), 
  'astrocyte' = c('GFAP', 'AQP4'),
  "Endo/Mural" = c("CLDN5", "CARMN", "ITIH5", "NOTCH3", "ATP10A", "MECOM", "EBF1", 
                   "AC092957.1", "ITGA1", "VWF"),
  "Choroid Plexus" = c("klotho", "CLIC6", "OATP14", "Ezrin")
)

# changing rownames for gene annotation purposes 
rownames(sce) <- rowData(sce)$Symbol

####### PLOTTING ###############################################################
# creating plot dir
plot_dir <- here("plots", "05_explore_sce", "03_Plot_Investigate_snAnno", "01_snAnno_orig_Violin")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

# creating name of pdf
snAnnoCustom <- here(plot_dir, "snAnno_custom_markers_violin_plots.pdf")

# plotting by snAnno in sce because snAnno is my original combined cluster resolution for sn annotations
my_plotMarkers(sce, marker_list = markers.custom, assay = "logcounts", 
               cat = "snAnno", fill_colors = NULL, pdf_fn = snAnnoCustom)

# tell Louise that sometimes code can get hooked when missing markers. Fix is to cancel function and then dev.off()
