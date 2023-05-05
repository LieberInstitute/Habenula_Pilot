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
library(tidyverse)

# loading final sce object 
load(file = here("processed-data", "99_paper_figs", "sce_objects", 
                 "official_final_sce.RDATA"), verbose = TRUE)
# sce

# loading final rse object (didn't change the location for this one)
load(here("processed-data", "02_bulk_qc", "count_data_bukola", 
          "rse_gene_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata"), verbose = TRUE)
# rse_gene

# creating plotting directory
plot_dir <- here("plots", "99_paper_figs", "10_bulk_Deconvo", "OPC_noisy")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

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
                        pdf_fn = here(plot_dir, "Top_10_Markers_OPC_noisy.pdf")
)

## copied directly from run_Bisque.R file in the bulk deconvo folder 
## creating marker_list of top 25 genes
marker_genes <- marker_stats |>
  filter(rank_ratio <= 25, gene %in% rownames(rse_gene)) |>
  pull(gene)

length(marker_genes)
# [1] 170

##### Running BISQUE ###########################################################
exp_set_bulk <- Biobase::ExpressionSet(assayData = assays(rse_gene[marker_genes,])$counts,
                                       phenoData=AnnotatedDataFrame(
                                         as.data.frame(colData(rse_gene))[c("BrNum")]))

exp_set_sce <- Biobase::ExpressionSet(assayData = as.matrix(assays(sym_sce[marker_genes,])$counts),
                                      phenoData=AnnotatedDataFrame(
                                        as.data.frame(colData(sym_sce))[,c("bulkTypeSepHb","Sample")]))

# checking for nuclei with 0 marker expression
zero_cell_filter <- colSums(exprs(exp_set_sce)) != 0
message("Exclude ", sum(!zero_cell_filter), " cells")
# Exclude 0 cells [yayyyyy!!!!]

est_prop <- ReferenceBasedDecomposition(bulk.eset = exp_set_bulk,
                                        sc.eset = exp_set_sce,
                                        cell.types = "bulkTypeSepHb",
                                        subject.names = "Sample",
                                        use.overlap = FALSE)

#### grabbed from my exploreBisque.R file in the bulk deconvo folder
# custom color scheme
color_bulk_clusters <- 
  c( "Oligo" = c("#A9A9A9"), # dark grey
     "OPC"= c("#7393B3"), # blue grey
     "Microglia" = c("#E5E4E2"), # platinum
     "Astrocyte" = c("#36454F"), # ash grey
     "Endo" = c("#848884"), # smoke
     "Inhib.Thal" = c('#2AAA8A'), # jungle green
     "Excit.Thal" = c("#478778"), # lincoln green
     "LHb" = c("#DE3163"), # cerise
     "MHb" = c("#00FFFF") # aqua
  )

# grabbing relevant phenotype info for bulk data
pd <- colData(rse_gene) |>
  as.data.frame() |>
  select(Sample = RNum, BrNum, PrimaryDx)

est_prop$bulk.props <- t(est_prop$bulk.props)
head(est_prop$bulk.props)

prop_long <- est_prop$bulk.props |>
  as.data.frame() |>
  rownames_to_column("Sample") |>
  tidyr::pivot_longer(!Sample, names_to = "cellType", values_to = "prop") |>
  left_join(pd) |>
  mutate(factor_CT = factor(cellType, levels = 
                              c("Astrocyte", "Endo", "Microglia", "Oligo", "OPC",
                                "Inhib.Thal", "Excit.Thal" , "MHb", "LHb")) ) |>
  arrange(factor_CT)

sum_Prop <- prop_long |>
  filter(cellType %in% c("LHb", "MHb")) |>
  group_by(BrNum) |>
  summarize(Hb_sum = sum(prop)) |>
  mutate(Br_Order = fct_reorder(BrNum, Hb_sum)) |>
  arrange(Br_Order)

prop_long <- left_join(prop_long, sum_Prop) |>
  arrange(Br_Order) 


## create composition bar plots
pdf(here(plot_dir, "bulk_Deconvo_Composition_OPC_noisy.pdf"), width = 21, height = 12)
plot_composition_bar(prop_long = prop_long, sample_col = "Br_Order",
                     x_col = "Br_Order", ct_col = "factor_CT") + 
  scale_fill_manual(values = color_bulk_clusters) +
  ggtitle("OPC_noisy")
dev.off()

# 