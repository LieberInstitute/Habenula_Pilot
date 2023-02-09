## 2/9/23 - Bukola Ajanaku
# Annotating trails for my three top clustering methods.
# qrsh -l mem_free=20G,h_vmem=20G

library("dendextend")
library("dynamicTreeCut")
library("SingleCellExperiment")
library("batchelor")
library("scater")
library("scran")
library("uwot")
library("DropletUtils")
library("jaffelab")
library("Rtsne")
library("here")
library("utils")
library("sessioninfo")
library("scuttle")

# loading sce object with clustered harmonized data
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
     "sce_mid_clustering.Rdata"))

####### NOTE: ANNOTATING BY TWO TOP CLUSTERING METHODS #########################
# 1) "k_20_Erik" (Walk Trap 20: 22 Groups) Has Rand value (against regular ctErik) of 0.699.
# 2) "k_10_Erik" (WalkTrap 10: 36 Groups) Has Rand value (against regular ctErik) of 0.655.

# sourcing code from DLPFC Project (by Louise Huuki) for plotting gene marker expression
# good to source because renders images by fill rather than taking up memory by plotting
  # generates expression plot 
source(here("code", "04_snRNA-seq", "sourcing", "custom_plotExpression.R"))

  # actually runs and plots markers 
source(here("code", "04_snRNA-seq", "sourcing", "my_plotMarkers.R"))

####### Marker gene list #######################################################
# Markers from Josh were either redudant or not specific enough. Using Erik's 
# for now.
markers.custom = list(
  'neuron' = c('SYT1'),# 'SNAP25', 'GRIN1','MAP2'),
  'excitatory_neuron' = c('SLC17A6'),# 'SLC17A7', 'SLC17A8'),
  'inhibitory_neuron' = c('GAD1', 'GAD2'), #'SLC32A1'),
  'mediodorsal thalamus'= c('EPHA4','PDYN', 'LYPD6B', 'LYPD6', 'S1PR1', 'GBX2', 'RAMP3', 'COX6A2', 'SLITRK6', 'DGAT2'),
  'Hb neuron specific'= c('POU2F2','POU4F1','GPR151','CALB2'),#,'GPR151','POU4F1','STMN2','CALB2','NR4A2','VAV2','LPAR1'),
  'MHB neuron specific' = c('TAC1','CHAT','CHRNB4'),#'TAC3','SLC17A7'
  'LHB neuron specific' = c('HTR2C','MMRN1'),#'RFTN1'
  'oligodendrocyte' = c('MOBP'),# 'MBP', 'PLP1'),
  'oligodendrocyte_precursor' = c('PDGFRA'),# 'VCAN', 'CSPG4', 'GPR17'),
  'microglia' = c('C3'),# 'CSF1R', 'C3'),
  'astrocyte' = c('GFAP')#,# 'AQP4'),
)

# Prepping sce object to plot by gene expression ###############################
rownames(sce) <- rowData(sce)$Symbol
# adding logcounts 
sce <- logNormCounts(sce)

###### Plotting gene expression for walktrap method 20 (22 groups) #############
pdf20 <- here("plots", "04_snRNA-seq", "07_Gene_Marking", "wt20_annotations.pdf")

my_plotMarkers(sce, marker_list = markers.custom, assay = "logcounts", 
               cat = "k_20_Erik", fill_colors = NULL, pdf_fn = pdf20)


###### Plotting gene expression for walktrap method 10 (36 groups) #############
pdf10 <- here("plots", "04_snRNA-seq", "07_Gene_Marking", "wt10_annotations.pdf")

my_plotMarkers(sce, marker_list = markers.custom, assay = "logcounts", 
               cat = "k_10_Erik", fill_colors = NULL, pdf_fn = pdf10)

### Saving #####################################################################
save(sce, file = here("processed-data", "04_snRNA-seq", "sce_objects", 
                      "sce_clustered_predrop.Rdata"))

# sgejobs::job_single('07_Gene_Marking.R', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript get_droplet_scores.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
