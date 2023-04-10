## 2/9/23 - Bukola Ajanaku
# Using walkTrap 10, I'll be comparing DRD5 and CHAT expression!
# qrsh -l mem_free=20G,h_vmem=20G

library("SingleCellExperiment")
library("here")
library("ggplot2")
library("DeconvoBuddies")

# getting sce object
load(here("processed-data", "99_paper_figs",  "sce_objects", "paper_worthy_sce.RDATA"),
     verbose = TRUE)

# creating plot_dir
plot_dir <- here("plots", "99_paper_figs", "99_Kristen_DRD5")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

# checking rownames of sce for symbol
head(rownames(sce))

# sourcing code for my_plotMarkers
source(here("code", "04_snRNA-seq", "sourcing", "custom_plotExpression.R"))
source(here("code", "04_snRNA-seq", "sourcing", "my_plotMarkers.R"))

# marker.list
markers.custom <- list(
  "DRD5_vs_CHAT" <- c("DRD5", "CHAT")
)

# plotting
pdf(here(plot_dir, "DRD5vsCHAT.pdf"))
  my_plotMarkers(sce, marker_list = markers.custom, assay = "logcounts", 
               cat = "snAnno2", fill_colors = NULL)
dev.off()