## Based on
## https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/master/10x_all-FACS-n10_2021rev_step01_processing-QC_MNT.R

library("SingleCellExperiment")
library("DropletUtils")
# library("BiocParallel")
library("scuttle")
library("tidyverse")
library("here")
library("sessioninfo")

## get user supplied sample + lower cutoff
args <- commandArgs(trailingOnly = TRUE)
user_sample <- args[[1]]
user_lower <- as.integer(args[[2]])


Br<-unlist(str_split(user_sample, "/"))[9]
#### Load & Subset raw data ####
load(here("processed-data","09_snRNA-seq_re-processed","20220601_human_hb_processing.rda"))

stopifnot(user_sample %in% sce.all.hb$Sample)

message("Running Sample: ", Br)

sce.all.hb <- sce.all.hb[, sce.all.hb$Sample == user_sample]
message("ncol:", ncol(sce.all.hb))

#### Run barcodeRanks to find knee ####

bcRanks <- barcodeRanks(sce.all.hb, fit.bounds = c(10, 1e3))

knee_lower <- metadata(bcRanks)$knee + 100
message(
    "'Second knee point' = ", metadata(bcRanks)$knee, "\n",
    "knee_lower = ", knee_lower, "\n",
    "user defined lower = ", user_lower
)

#### Run emptyDrops w/ knee + 100 ####
set.seed(100)
message("Starting emptyDrops")
Sys.time()
e.out <- DropletUtils::emptyDrops(
    sce.all.hb,
    niters = 30000,
    lower = user_lower
    # ,
    # BPPARAM = BiocParallel::MulticoreParam(4)
)
message("Done - saving data")
Sys.time()

save(e.out, file = here("processed-data","09_snRNA-seq_re-processed", "droplet_scores_troubleshoot", paste0("droplet_scores_", Br, ".Rdata")))

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
    geom_hline(yintercept = user_lower, linetype = "dashed", color = "red") +
    annotate("text", x = 10, y = user_lower, label = "User defined 'lower'", vjust = -0.5, color = "red") +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    labs(
        x = "Barcode Rank",
        y = "Total UMIs",
        title = paste("Sample", Br),
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

ggsave(droplet_elbow_plot, filename = here("plots","09_snRNA-seq_re-processed", "droplet_knee_plots", paste0("droplet_qc_", Br, ".png")))


# sgejobs::job_single('get_droplet_scores_troubleshoot', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript get_droplet_scores_troubleshoot.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
