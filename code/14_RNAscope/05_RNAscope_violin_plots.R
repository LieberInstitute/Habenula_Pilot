
library("SingleCellExperiment")
library("here")
library("ggplot2")
library("DeconvoBuddies")

# loading sce object
load(here("processed-data",  "sce_objects", "official_final_sce.RDATA"),
     verbose = TRUE)
# sce

# creating plot directory
plot_dir <- here("plots", "14_RNAscope", "05_RNAscope_violin_plots")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

# sourcing official color palette
source(file = here("code", "99_paper_figs", "source_colors.R"))
# bulk_colors and sn_colors

# using Louise's package (yayyyy!!) to make the violin plots easy peasy
rownames(sce) <- rowData(sce)$Symbol

LHb1 <- plot_gene_express(sce = sce,
                          genes = c("POU4F1", "SEMA3D", "TLE2", "ONECUT2"),
                          cat = "final_Annotations",
                          color_pal = sn_colors,
                          title = "Expression in snRNA-seq",
                          ncol = 4)

ggsave(LHb1, filename = here(plot_dir, "RNAscope_violin_LHb1.pdf"),  height = 4, width = 9)


LHb2 <- plot_gene_express(sce = sce,
                          genes = c("POU4F1", "CRH", "MCOLN3","ESRP1"),
                          cat = "final_Annotations",
                          color_pal = sn_colors,
                          title = "Expression in snRNA-seq",
                          ncol = 4)

ggsave(LHb2, filename = here(plot_dir, "RNAscope_violin_LHb2.pdf"),  height = 4, width = 9)


MHb1 <- plot_gene_express(sce = sce,
                          genes = c("POU4F1","CHAT", "EBF3","CCK"),
                          cat = "final_Annotations",
                          color_pal = sn_colors,
                          title = "Expression in snRNA-seq",
                          ncol = 4)

ggsave(MHb1, filename = here(plot_dir, "RNAscope_violin_MHb1.pdf"),  height = 4, width = 9)

MHb2 <- plot_gene_express(sce = sce,
                          genes = c("CHRNB4","CHAT","CCK"),
                          cat = "final_Annotations",
                          color_pal = sn_colors,
                          title = "Expression in snRNA-seq",
                          ncol = 4)

ggsave(MHb2, filename = here(plot_dir, "RNAscope_violin_MHb2.pdf"),  height = 4, width = 7)

