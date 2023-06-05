## June 5, 2023 - Bukola Ajanaku
# Working on MAGMA trans-special analysis of our sce object against the 
# Wallace et al. 2019 paper mouse data.
# qrsh -l mem_free=50G,h_vmem=50G

library("SingleCellExperiment")
library("here")
library("sessioninfo")
library("jaffelab")

# loading our final sce object
load(file = here("processed-data", "99_paper_figs", "sce_objects", 
                 "official_final_sce.RDATA"), verbose = TRUE)
  # sce 

# checking cluster data 
table(sce$final_Annotations)

# loading Wallace et al. clean data
load(file = here("processed-data", "99_paper_figs", "MAGMA",
                 "Wallace_mouse_data.rds"), verbose = TRUE)
  # 

# grabbed from Matt's code 