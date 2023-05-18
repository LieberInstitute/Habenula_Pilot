## May 18, 2023 - Bukola Ajanaku
# Grabbing marker genes for sn-RNA seq resolution of clusters.
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

table(sce$final_Annotations)
    # Astrocyte       Endo Excit.Thal Inhib.Thal      LHb.1      LHb.2      LHb.3 
    # 538         38       1800       7612        201        266        134 
    # LHb.4      LHb.5      LHb.6      LHb.7      MHb.1      MHb.2      MHb.3 
    # 477         83         39       1014        152        540         18 
    # Microglia      Oligo        OPC 
    # 145       2178       1202 

# creating plotting directory
plot_dir <- here("plots", "99_paper_figs", "10c_snResolution_Top_Markers")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

# sourcing official color palette 
source(file = here("code", "99_paper_figs", "source_colors.R"))
# bulk_colors and sn_colors

# adding necessary labels 
sym_sce <- sce
rownames(sym_sce) <- rowData(sce)$Symbol

# pre-bisque
# Creating mean_ratios based on our specified annotations
ratios <- get_mean_ratio2(sym_sce,
                          cellType_col = "final_Annotations",
                          assay_name = "logcounts",
                          add_symbol = TRUE)

# Using the 1 vs All standard fold change for each gene x cell type
fc <- findMarkers_1vAll(sym_sce,
                        assay_name = "counts",
                        cellType_col = "final_Annotations",
                        add_symbol = FALSE,
                        mod = "~Sample",
                        verbose = TRUE
)

# combining the two to form marker_stats
marker_stats <- left_join(ratios, fc, by = c("gene", "cellType.target"))

######### PLOTTING TOP 10 GENE MARKERS #########################################
plot_marker_express_ALL(sym_sce,
                        marker_stats,
                        n_genes = 15,
                        rank_col = "rank_ratio",
                        anno_col = "anno_ratio",
                        cellType_col = "final_Annotations",
                        pdf_fn = here(plot_dir, "snResolution_top_15_MarkerGenes.pdf") 
                        
) 

########## PLOTTING HOCKEY STICKS ##############################################
# adding color group
marker_stats$Top25 <- "No"
marker_stats[which(marker_stats$rank_ratio <= 25), "Top25"] <- "Yes"

# plotting hockey sticks 
pdf(here(plot_dir, "snResolution_hockeysticks_top25.pdf"))
  pos = position_jitter(seed = 1)
  
  ggplot(marker_stats, aes(ratio, std.logFC)) +
  geom_jitter(size = 0.5, 
               aes(colour = Top25),
               position = pos) + 
    labs(x = "Mean Ratio") +
    guides(colour = guide_legend(title = "Top 25 Marker")) +
    geom_text(
      data = subset(marker_stats, Top25 == "Yes"),
      position = pos,
      aes(
        label = Symbol
      )
    ) 
  
dev.off()

########### EXPORTING TOP MARKER GENES #########################################
exp_Mark_Table <- marker_stats |>
  filter(rank_ratio <= 50)

# exporting (into plotting directory (i know)) as a csv table
write.xlsx(exp_Mark_Table, file = here(plot_dir, "snResolution_top50MarkerGenes.xlsx"),
           sheetName = "snRNA-seq Resolution Top 50 Marker Genes", append = FALSE)









# 