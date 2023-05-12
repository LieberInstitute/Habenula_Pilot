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
library(xlsx)
library(sessioninfo)

# loading final sce object 
load(file = here("processed-data", "99_paper_figs", "sce_objects", 
                 "official_final_sce.RDATA"), verbose = TRUE)
# sce
dim(sce)
  # [1] 33848 16437 (dropped OPC_noisy and Excit.Neuron ambig)

# loading final rse object (didn't change the location for this one)
load(here("processed-data", "02_bulk_qc", "count_data_bukola", 
          "rse_gene_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata"), verbose = TRUE)
# rse_gene

# creating plotting directory
plot_dir <- here("plots", "99_paper_figs", "10_bulk_Deconvo", "OPC_clean")
  # if(!dir.exists(plot_dir)){
  #   dir.create(plot_dir)
  # }

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

    ## grabbing the top 50 genes for export 
    # exp_Mark_Table <- marker_stats |>
    #   filter(rank_ratio <= 50)
    
    # exporting (into plotting directory (i know)) as a csv table
    # write.xlsx(exp_Mark_Table, file = here(plot_dir, "top50MarkerGenes.xlsx"),
    #            sheetName = "Top 50 Marker Genes Per Cell Type", append = FALSE)
    # 


######### PLOTTING TOP 10 GENE MARKERS #########################################
plot_marker_express_ALL(sym_sce,
                        marker_stats,
                        n_genes = 10,
                        rank_col = "rank_ratio",
                        anno_col = "anno_ratio",
                        cellType_col = "bulkTypeSepHb",
                        pdf_fn = here(plot_dir, "Top_10_Markers_OPC_clean.pdf")
)
########## PLOTTING HOCKEY STICKS ##############################################
# adding color group
marker_stats$Top25 <- "No"
marker_stats[which(marker_stats$rank_ratio <= 25), "Top25"] <- "Yes"

# plotting hockey sticks 
pdf(here(plot_dir, "hockeysticks_Official.pdf"))
ggplot(marker_stats, aes(ratio, std.logFC)) +
  geom_point(size = 0.5, aes(colour = Top25)) +  
  facet_wrap(~cellType.target, scales = "free") +
  labs(x = "Mean Ratio") +
  guides(colour = guide_legend(title = "Top 25 Marker"))
dev.off()


##### Running BISQUE ###########################################################
## creating marker_list of top 25 genes
marker_genes <- marker_stats |>
  filter(rank_ratio <= 25, gene %in% rownames(rse_gene)) |>
  pull(gene)

length(marker_genes)
# [1] 170

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
  arrange(Br_Order) |>
  mutate(prop_perc = prop * 100)

## create composition bar plots
pdf(here(plot_dir, "bulk_Deconvo_Composition_OPC_clean.pdf"), width = 21, height = 11)
comp_plot <- ggplot(prop_long, 
                    aes(x = Br_Order, y = prop_perc, fill = factor_CT)) +
  geom_col() +
  geom_text(
    data = subset(prop_long, prop_perc >= 1),
    aes(
      label = round(prop_perc, 1)
    ), 
    size = 4,
    position = position_stack(vjust = 0.4),
    angle = -90,
    fontface = "bold",
    colour = "black"
  ) +
  scale_fill_manual(values = alpha(color_bulk_clusters, 0.8)) +
  theme_bw() +
  theme(legend.position = "None", 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(size = 16, hjust = 0.5)) +
  ggtitle("Bulk Deconvolution") +
  ylab("Proportion")

comp_plot

dev.off()
  
  # pdf(here(plot_dir, "bulk_Deconvo_Composition_OPC_clean.pdf"), width = 21, height = 12)
  #   plot_composition_bar(prop_long = prop_long, sample_col = "Br_Order",
  #                        x_col = "Br_Order", ct_col = "factor_CT") + 
  #     scale_fill_manual(values = color_bulk_clusters) +
  #     ggtitle("Bulk Deconvolution") + 
  #     theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
  #     
  # dev.off()

  
  
  

#