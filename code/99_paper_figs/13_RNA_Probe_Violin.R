## June 15, 2023 - Bukola Ajanaku
# This is the exact code I used for the original progress report. Copy and pasted!
# Simply changed plotting directory!!!
# qrsh -l mem_free=30G,h_vmem=30G

# loading relevant libraries
library("SingleCellExperiment")
library("here")
library("ggplot2")
library("DeconvoBuddies")

# loading sce object 
load(here("processed-data", "04_snRNA-seq",  "sce_objects", "sce_final.Rdata"),
     verbose = TRUE)
# sce_final 

# creating plot directory
plot_dir <- here("plots", "99_paper_figs", "13_RNA_Probe_Violin")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

# sourcing official color palette 
source(file = here("code", "99_paper_figs", "source_colors.R"))
# bulk_colors and sn_colors

marker_list <- c(
  MHb = c("TAC3", "CCK", "CHRNB4", "BHLHE22", "CHAT"),
  LHb = c("ONCUT2", "MCOLN3", "SEMA3D", "HTR4"),
  LHb = c("ESRP1", "TLE2", "CRH", "HTR4")
)


# using Louise's package (yayyyy!!) to make the violin plots easy peasy
plot_marker_express_List(sce = sce_final, 
                         gene_list = marker_list,
                         cellType_col = "final_Annotations",
                         gene_name_col = "Symbol",
                         color_pal = sn_colors,
                         pdf_fn = here(plot_dir, "RNA_probe_Violin_Plot.pdf"))