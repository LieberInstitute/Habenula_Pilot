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
          "sce_post_clustering_with_logcounts.Rdata"))

####### NOTE: ANNOTATING BY TWO THREE CLUSTERING METHODS #######################
# Erik used Walktrap 50 and ended up with 14 groups.
# 1) Walktrap 50 (Rand: 0.661) has 14 groups.
# 2) Walktrap 10 (Rand: 0.654) has 37 groups.
# 3) Walktrap 20 (Rand: 0.632) has 23 groups.

####### SOURCING ###############################################################
# sourcing code from DLPFC Project (by Louise Huuki) 
# plots gene expression in a manner that renders images by fill rather than taking 
# up memory by plotting each point
  # generates expression plot 
source(here("code", "04_snRNA-seq", "sourcing", "custom_plotExpression.R"))

  # actually runs and plots markers 
source(here("code", "04_snRNA-seq", "sourcing", "my_plotMarkers.R"))

####### Marker gene list #######################################################
# [MUST HAVE AT LEAST TWO GENES PER MARKER FOR FUNCTION TO WORK CORRECTLY]
  # Erik's markers
    # markers.custom = list(
    #   'neuron' = c('SYT1', 'SNAP25'), #'SNAP25', 'GRIN1','MAP2'),
    #   'excitatory_neuron' = c('SLC17A6', 'SLC17A7'), # 'SLC17A8'),
    #   'inhibitory_neuron' = c('GAD1', 'GAD2'), #'SLC32A1'),
    #   'mediodorsal thalamus'= c('EPHA4','PDYN', 'LYPD6B', 'LYPD6', 'S1PR1', 'GBX2', 'RAMP3', 'COX6A2', 'SLITRK6', 'DGAT2'),
    #   'Hb neuron specific'= c('POU2F2','POU4F1','GPR151','CALB2'),#,'GPR151','POU4F1','STMN2','CALB2','NR4A2','VAV2','LPAR1'),
    #   'MHB neuron specific' = c('TAC1','CHAT','CHRNB4'),#'TAC3','SLC17A7'
    #   'LHB neuron specific' = c('HTR2C','MMRN1'),#'RFTN1'
    #   'oligodendrocyte' = c('MOBP', 'MBP'), # 'PLP1'),
    #   'oligodendrocyte_precursor' = c('PDGFRA', 'VCAN'), # 'CSPG4', 'GPR17'),
    #   'microglia' = c('C3', 'CSF1R'), #'C3'),
    #   'astrocyte' = c('GFAP', 'AQP4')
    # )

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

#### PREPPING sce object to plot by gene expression ############################
# adding logcounts 
# sce <- logNormCounts(sce)
# 
# message("Start - saving data")
# Sys.time()
# 
# # saving before messing with row names for annotation purposes
# save(sce, file = here("processed-data", "04_snRNA-seq", "sce_objects",
#                         "sce_post_clustering_with_logcounts.Rdata"))
# message("Done - saving data")
# Sys.time()

# changing rownames for gene annotation purposes 
rownames(sce) <- rowData(sce)$Symbol


###### Plotting gene expression for walktrap method 10 (37 groups) #############
message("Start - Annotating wt10")
  Sys.time()
  
pdf10 <- here("plots", "04_snRNA-seq", "07_Gene_Marking", "wT_10_annotations_more_gran1.pdf")

my_plotMarkers(sce, marker_list = markers.custom, assay = "logcounts", 
               cat = "wT_10_Erik", fill_colors = NULL, pdf_fn = pdf10)

message("End - Annotating wt10")
Sys.time()

###### Plotting gene expression for walktrap method 20 (23 groups) #############
message("Start - Annotating wt20")
Sys.time()

pdf20 <- here("plots", "04_snRNA-seq", "07_Gene_Marking", "wT_20_annotations_more_gran1.pdf")

my_plotMarkers(sce, marker_list = markers.custom, assay = "logcounts", 
               cat = "wT_20_Erik", fill_colors = NULL, pdf_fn = pdf20)

message("End - Annotating wt20")
Sys.time()

###### Plotting gene expression for walktrap method 50 (14 groups) #############
message("Start - Annotating wt50")
  Sys.time()
  
  pdf50 <- here("plots", "04_snRNA-seq", "07_Gene_Marking", "wT_50_annotations_more_gran_1.pdf")
  
  my_plotMarkers(sce, marker_list = markers.custom, assay = "logcounts", 
                 cat = "wT_50_Erik", fill_colors = NULL, pdf_fn = pdf50)
  
  message("End - Annotating wt50")
Sys.time()



# sgejobs::job_single('07_Gene_Marking', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript 07_Gene_Marking.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
