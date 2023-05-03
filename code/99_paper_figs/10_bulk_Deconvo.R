## May 3, 2023 - Bukola Ajanaku
# Running Bulk Deconvolution Plots for Top 10 gene markers! 
# qrsh -l mem_free=50G,h_vmem=50G

library(SummarizedExperiment)
library(here)
library(DeconvoBuddies)
library(SingleCellExperiment)
library(jaffelab)
library(dplyr)
library(BisqueRNA)
library(ggplot2)

# loading final sce object 
load(file = here("processed-data", "99_paper_figs", "sce_objects", 
                 "official_final_sce.RDATA"), verbose = TRUE)
# sce

# loading final rse object (didn't change the location for this one)
load(here("processed-data", "02_bulk_qc", "count_data_bukola", 
          "rse_gene_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata"), verbose = TRUE)
# rse_gene

# creating plotting directory
plot_dir <- here("plots", "99_paper_figs", "10_bulk_Deconvo", "OPC_clean")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

### CLEANING OPC ###############################################################
dim(sce)
  # [1] 33848 17031

sce <- sce[, which(sce$OPC_clean == "Yes")]

dim(sce)
  # [1] 33848 16437

###### Adding bulk collapsed annotations to sce object #########################
sce$bulkTypeSepHb <- sce$final_Annotations
# making separated Hb (2)
sce$bulkTypeSepHb[sce$bulkTypeSepHb %in% grep("^LHb\\.", unique(sce$bulkTypeSepHb), value = TRUE)] <-  'LHb'
sce$bulkTypeSepHb[sce$bulkTypeSepHb %in% grep("^MHb\\.", unique(sce$bulkTypeSepHb), value = TRUE)] <-  'MHb'

# check!
table(sce$bulkTypeSepHb)
    # Astrocyte   Endo Excit.Thal Inhib.Thal        LHb        MHb  Microglia 
    # 538         38       1800       7612       2214        710        145 
    # Oligo        OPC 
    # 2178       1202 

###### Adding necessary symbols ################################################
sym_sce <- sce
rownames(sym_sce) <- rowData(sce)$Symbol

rownames(rse_gene) <- rowData(rse_gene)$Symbol

######## Pre-Bisque ############################################################
## remember, this is the broad analyses meaning that these annotations are solely
# for bulk deconvo

# Creating mean_ratios based on our specified annotations
ratios <- get_mean_ratio2(sym_sce,
                          cellType_col = "bulkTypeSepHb",
                          assay_name = "logcounts",
                          add_symbol = TRUE)

# Using the 1 vs All standard fold change for each gene x cell type
fc <- findMarkers_1vAll(sym_sce,
                        assay_name = "counts",
                        cellType_col = "bulkTypeSepHb",
                        add_symbol = FALSE,
                        mod = "~Sample",
                        verbose = TRUE
)

# combining the two to form marker_stats
marker_stats <- left_join(ratios, fc, by = c("gene", "cellType.target"))

# Random color scheme [NEED TO ESTABLISH MY OWN FOR THIS STEP]
cell_types <- levels(sym_sce$cellType)
# cell_colors <- create_cell_colors(cell_types = cell_types, pallet = "classic", split = "\\.", preview = TRUE)
# cell_colors errorr

# printing top 10 markers 
plot_marker_express_ALL(sym_sce,
                        marker_stats,
                        n_genes = 10,
                        rank_col = "rank_ratio",
                        anno_col = "anno_ratio",
                        cellType_col = "bulkTypeSepHb",
                        pdf_fn = here(plot_dir, "Top_10_Markers_OPC_clean.pdf")
)





# 