## April 11, 2023 - Bukola Ajanaku
# Working on bulkRNA-seq deconvolution using Bisque.
# qrsh -l mem_free=50G,h_vmem=50G

library(SummarizedExperiment)
library(here)
library(DeconvoBuddies)
library(SingleCellExperiment)
library(jaffelab)
library(dplyr)
library(BisqueRNA)
library(ggplot2)

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
    # 2178       1202


# adding symbols
sym_sce <- sce
rownames(sym_sce) <- rowData(sce)$Symbol

# check 
dim(sym_sce)
  # [1] 33848 16437


# adding symbols 
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
                      pdf_fn = here(plot_dir, "Top_10_Broad_snRNA_Annotation_Markers.pdf")
)

# OPC is pretty messy. Only taking top markers 1-6! 
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
##### NOTE: DROPPING 0PC Samples ###############################################
# OPC cluster is very dirty and upon inspection of previous heatmaps, three 
# Samples are contributing to the noise in this cluster and will be dropped.
# They are as follows: "Br5555", "Br1204", and "Br1092".

# I haven't edited the previous copies of this sce. Will need to go back and drop 
# where I dropped the Hb cluster!

# adding color group
marker_stats$Top25 <- "No"
marker_stats[which(marker_stats$rank_ratio <= 25), "Top25"] <- "Yes"

# plotting
# this is the hockey stick plot for the split up snAnno
pdf(here(plot_dir, "hockystick_cleanedOPC.pdf"))
  ggplot(marker_stats, aes(ratio, std.logFC)) +
    geom_point(size = 0.5, aes(colour = Top25)) +  
    facet_wrap(~cellType.target, scales = "free") 
dev.off()


#### Saving
save(marker_stats, file = here(save_dir, "marker_stats_top_25_genes.csv"))

save(est_prop, file = here(save_dir, "est_prop_split_Hb_annotations.RDATA"))


# 