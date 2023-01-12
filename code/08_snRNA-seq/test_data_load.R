
library("here")
library("SingleCellExperiment")
library("DropletUtils")

#### load the whole data

message(Sys.time(), "- load full sce...")

load_start <- Sys.time()
load(here("processed-data","09_snRNA-seq_re-processed","20220601_human_hb_processing.rda"), verbose = TRUE)
load_end <- Sys.time()

message("Time to load sce")
load_end - load_start

# ncol(sce)

#### read10xCounts
samples <- list.files(here("processed-data","07_cellranger"))

message(Sys.time() - "Use read10xCounts")

read_10x_start <- Sys.time()
# Br5558/outs/raw_feature_bc_matrix
sce_10x <- read10xCounts(here("processed-data","07_cellranger", "Br5558", "outs", "raw_feature_bc_matrix"))
read_10x_end <- Sys.time()

message("Time to load sce")
read10x_end - read10x_start

# sgejobs::job_single('test_data_load', queue = 'bluejay', create_shell = TRUE, memory = '50G', command = "Rscript test_data_load.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
