## March 8 2023 - Bukola Ajanaku
# This folder is solely used to investigate the (19) single nucleus annotation clusters
# to validate probe their identities and investigate if any further collapses are necessary.

# Notes 3/9/23:
# 1) Plotting violin plots for combined clusters against custom marks. Noticed no real endo cells in Astrocytes
 # but noticed them in LHb.6.
 # 2) Completed heatmap that psuedobulked by snAnno classificatios and heatmapped it. [May need to update custom 
 #                                                                                     markers based on list from Maynard.]
 #    To Do: a) work on custom markers list
#           a) merge MHb 2&3 and re-run mean ratio to see how different everything looks in the violin format.
 #            ii) run heatmap expression on pseudo_bulked new snAnno.
 #            iii) re-run original violin plots on new snAnno.
 #           b) Work on hockey stick plots for :
 #            i) original snAnno
 #            ii) new snAnnno

# Official Habenula Markers List:
# 
new_markers.custom <- c(
'Neuron' = c('SYT1', 'SNAP25', "SYT4", "SYP"), 
"Non-Neuronal Subtype" = c("TnF", "Kir4.1", "Kcnj10"),
'excitatory_neuron' = c('SLC17A6', 'SLC17A7'), 
'inhibitory_neuron' = c('GAD1', 'GAD2'), 
"Habenula_neurons" = c("Ccbp2","Cd63", "Htr5b", "Kcnh8", "Kctd8", "Lrrc55", "Mapk4", "Neurod1", 
              "Pixnc1", "Scube1", "Sstr4", "Tacr1", "Sstr2", "Irx2"),
"LHB_neuron_specific" = c("Pcdh10", "Gabra1", "Syn2", "Gap43", "Htr2c", "Adcyap1r1", 
                     "Chrm3", "Vgf", "Gpr151", "Sst", "Slc17ar"),
"MHb_neuron_specific" = c("Tac3", "Spon1", "Sema3d", "Calb1",'Htr5b', "Cnr1", "Gpr4", "Chrnb3", "B4",
                   "Slc17a7"),
'mediodorsal thalamus'= c('EPHA4','PDYN', 'LYPD6B', 'LYPD6', 'S1PR1', 'GBX2', 
                          'RAMP3', 'COX6A2', 'SLITRK6', 'DGAT2', "ADARB2"),
"Endo/Mural" = c("CLDN5", "CARMN", "ITIH5", "NOTCH3", "ATP10A", "MECOM", "EBF1", 
                 "AC092957.1", "ITGA1", "VWF"),
"Choroid Plexus" = c("klotho", "CLIC6", "OATP14", "Ezrin"),
'oligodendrocyte' = c('MOBP', 'MBP', "Cx3cr1"),
'oligodendrocyte_precursor' = c('PDGFRA', 'VCAN'), 
'microglia' = c('C3', 'CSF1R'), 
'astrocyte' = c('GFAP', 'AQP4')
)

extra_markers.custom <- c(
"Macrophages" = c("Mfc1", "Mfc1"),
"Fibroblasts" = c("Pdgfra", "Col3a1"),
"Ependymal" = c("Slc6a11", "Hdc"),
"Pericytes" = c("Abcc9", "Pdgfrb"),
"Polydendro" =  c("Gpr17", "Olig1", 'Gap43', 'Pdgfra'), 

)


eLife_markers.custom <- c(
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








# 
