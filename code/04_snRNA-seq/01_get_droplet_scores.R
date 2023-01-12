# January 10, 2023 - Bukola Ajanaku
# Collaborated with Louise Huuki 
# Based on
# https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/master/10x_all-FACS-n10_2021rev_step01_processing-QC_MNT.R
# list.files("/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/07_cellranger")
# qrsh -l mem_free=75G,h_vmem=75G

library("SingleCellExperiment")
library("DropletUtils")
# library("BiocParallel")
library("scuttle")
library("tidyverse")
library("here")
library("sessioninfo")

## get sample i
args <- commandArgs(trailingOnly = TRUE)
sample <- args[[1]]
sample_path <- here("processed-data", "07_cellranger", sample, "outs", "raw_feature_bc_matrix") 
stopifnot(file.exists(sample_path))

message(Sys.time(), " Reading data from ", sample_path)

#### Load & Subset raw data ####
sce <- read10xCounts(sample_path, col.names=TRUE) 

message("ncol:", ncol(sce))

#### Run barcodeRanks to find knee ####
message(Sys.time(), "Running barcode ranks.")

bcRanks <- barcodeRanks(sce, fit.bounds = c(10, 1e3))

knee_lower <- metadata(bcRanks)$knee + 100
message(
  "'Second knee point' = ", metadata(bcRanks)$knee, "\n",
  "knee_lower =", knee_lower
)

#### Run emptyDrops w/ knee + 100 ####
set.seed(100)
message(Sys.time(), "Starting emptyDrops")
e.out <- DropletUtils::emptyDrops(
  sce,
  niters = 30000,
  lower = knee_lower
  # ,
  # BPPARAM = BiocParallel::MulticoreParam(4)
)
message(Sys.time(), "Done - saving data")

save(e.out, file = here("processed-data", "04_snRNA-seq_Bukola", "01_get_droplet_scores", paste0("droplet_scores_", sample, ".Rdata")))

#### QC Plots ####
message("QC check")
FDR_cutoff <- 0.001
addmargins(table(Signif = e.out$FDR <= FDR_cutoff, Limited = e.out$Limited, useNA = "ifany"))

n_cell_anno <- paste("Non-empty:", sum(e.out$FDR < FDR_cutoff, na.rm = TRUE))
message(n_cell_anno)

my_theme <- theme_bw() +
  theme(text = element_text(size = 15))

droplet_elbow_plot <- as.data.frame(bcRanks) %>%
  add_column(FDR = e.out$FDR) %>%
  ggplot(aes(x = rank, y = total, color = FDR < FDR_cutoff)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_hline(yintercept = metadata(bcRanks)$knee, linetype = "dotted", color = "gray") +
  annotate("text", x = 10, y = metadata(bcRanks)$knee, label = "Second Knee", vjust = -1, color = "gray") +
  geom_hline(yintercept = knee_lower, linetype = "dashed") +
  annotate("text", x = 10, y = knee_lower, label = "Knee est 'lower'", vjust = -0.5) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  labs(
    x = "Barcode Rank",
    y = "Total UMIs",
    title = paste("Sample", sample),
    subtitle = n_cell_anno,
    color = paste("FDR <", FDR_cutoff)
  ) +
  my_theme +
  theme(legend.position = "bottom")

# droplet_scatter_plot <- as.data.frame(e) %>%
#   ggplot(aes(x = Total, y = -LogProb, color = FDR < FDR_cutoff)) +
#   geom_point(alpha = 0.5, size = 1) +
#   labs(x = "Total UMIs", y = "-Log Probability",
#        color = paste("FDR <", FDR_cutoff)) +
#   my_theme+
#   theme(legend.position = "bottom")
# # print(droplet_elbow_plot/droplet_scatter_plot)
# ggsave(droplet_elbow_plot/droplet_scatter_plot, filename = here("plots","03_build_sce", "droplet_qc_png",paste0("droplet_qc_",sample,".png")))

ggsave(droplet_elbow_plot, filename = here("plots", "04_snRNA-seq_Bukola", "01_get_droplet_scores", paste0("droplet_qc_", sample, ".png")))


# sgejobs::job_single('01_get_droplet_scores', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript 01_get_droplet_scores.R Br1092")
# sample_list <- list.files("/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/07_cellranger")
# sgejobs::job_loop(loops = list(sample = sample_list), name = "01_get_droplet_scores_loop", 
#          create_shell = TRUE, queue = "bluejay", memory = "50G")



## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()