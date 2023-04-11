## April 11, 2023 - Bukola Ajanaku
# Working on bulkRNA-seq deconvolution 
# qrsh -l mem_free=20G,h_vmem=20G

library(SummarizedExperiment)
library(here)
library(DeconvoBuddies)
library(SingleCellExperiment)
library(jaffelab)
library(dplyr)
library(BisqueRNA)

# loading sce object (sn data) post drop of Hb cluster! 
load(here("processed-data", "06_deconvolution", "sce_objects", "sce_first_bulkTypes.RDATA"),
     verbose = TRUE)
      # sce 
# loading cleaned rse object (bulk data)
load(here("processed-data", "02_bulk_qc", "count_data_bukola", 
          "rse_gene_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata"), verbose = TRUE)
      # rse_gene

# creating plotting dir 
plot_dir <- here("plots", "06_deconvolution", "03_bulk_Composition")
  if(!dir.exists(plot_dir)){
    dir.create(plot_dir)
  }

# creating directory for processed data
save_dir <- here("processed-data", "06_deconvolution", "run_Bisque")
if(!dir.exists(save_dir)){
  dir.create(save_dir)
}

###### CHECKING DATA ###########################################################
table(sce$bulkTypeSepHb)
    # Astrocyte   Endo Excit.Thal   Inhib.Thal   LHb        MHb  Microglia 
    # 538         38       1800       7612       2214        710        145 
    # Oligo       OPC 
    # 2178       1796

# fixing sample issue
sce$Sample <- sce$RealSample
sce$FakeSample <- NULL
sce$RealSample <- NULL

# changing rownames of sce to actual gene symbols and saving as sce_symbol
sce_symbol <- sce
rownames(sce_symbol) <- rowData(sce)$Symbol

# changing rownames of rse_gene to actual gene symbols
rownames(rse_gene) <- rowData(rse_gene)$Symbol

######## Pre-Bisque ############################################################
## remember, this is the broad analyses meaning that these annotations are solely
# for bulk deconvo

# Creating mean_ratios based on our specified annotations
ratios <- get_mean_ratio2(sce_symbol,
                          cellType_col = "bulkTypeSepHb",
                          assay_name = "logcounts",
                          add_symbol = TRUE)

# Using the 1 vs All standard fold change for each gene x cell type
fc <- findMarkers_1vAll(sce_symbol,
                        assay_name = "counts",
                        cellType_col = "bulkTypeSepHb",
                        add_symbol = FALSE,
                        mod = "~Sample",
                        verbose = TRUE
                       )

# combining the two to form marker_stats
marker_stats <- left_join(ratios, fc, by = c("gene", "cellType.target"))

# Random color scheme [NEED TO ESTABLISH MY OWN FOR THIS STEP]
cell_types <- levels(sce_symbol$cellType)
# cell_colors <- create_cell_colors(cell_types = cell_types, pallet = "classic", split = "\\.", preview = TRUE)
  # cell_colors errorr

# printing top 10 markers 
plot_marker_express_ALL(sce_symbol,
                      marker_stats,
                      n_genes = 10,
                      rank_col = "rank_ratio",
                      anno_col = "anno_ratio",
                      cellType_col = "bulkTypeSepHb",
                      pdf_fn = here(plot_dir, "Top_10_Broad_snRNA_Annotation_Markers.pdf")
)

# creating marker_list of top 25 genes
marker_genes <- marker_stats |>
  filter(rank_ratio <= 25, gene %in% rownames(rse_gene)) |>
  pull(gene)

length(marker_genes)
# [1] 172

##### Running BISQUE ###########################################################
exp_set_bulk <- Biobase::ExpressionSet(assayData = assays(rse_gene[marker_genes,])$counts,
                                       phenoData=AnnotatedDataFrame(
                                         as.data.frame(colData(rse_gene))[c("BrNum")]))

exp_set_sce <- Biobase::ExpressionSet(assayData = as.matrix(assays(sce_symbol[marker_genes,])$counts),
                                      phenoData=AnnotatedDataFrame(
                                        as.data.frame(colData(sce_symbol))[,c("bulkTypeSepHb","Sample")]))

# checking for nuclei with 0 marker expression
zero_cell_filter <- colSums(exprs(exp_set_sce)) != 0
message("Exclude ", sum(!zero_cell_filter), " cells")
# Exclude 0 cells [yayyyyy!!!!]

est_prop <- ReferenceBasedDecomposition(bulk.eset = exp_set_bulk,
                                        sc.eset = exp_set_sce,
                                        cell.types = "bulkTypeSepHb",
                                        subject.names = "Sample",
                                        use.overlap = FALSE)

#### Saving
save(marker_stats, file = here(save_dir, "marker_stats_top_25_genes.csv"))

save(est_prop, file = here(save_dir, "est_prop_split_Hb_annotations.RDATA"))








# 