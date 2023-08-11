library("here")
library("data.table")
library("VennDiagram")
library("dplyr")
library("sessioninfo")

out_plot <- here("plots", "10_DEA", "08_Comparisions")
if (!dir.exists(out_plot)) dir.create(out_plot)
out_data <- here("processed-data", "10_DEA", "08_Comparisions")
if (!dir.exists(out_data)) dir.create(out_data)



############################### Load Sig DE data ##############################

message(Sys.time(), " - loading objects")

sig_genes <- fread(
    here("processed-data/10_DEA/04_DEA/DEA_Sig-gene_FDR1_qc-totAGene-qSVs-Hb-Thal.tsv"),
    sep = "\t",
    data.table = FALSE
)

sig_jx <- fread(
    here("processed-data/10_DEA/04_DEA/DEA_Sig-jx_FDR1_qc-totAGene-qSVs-Hb-Thal.tsv"),
    sep = "\t",
    data.table = FALSE
)


sig_exons <- fread(
    here("processed-data/10_DEA/04_DEA/DEA_Sig-exon_FDR1_qc-totAGene-qSVs-Hb-Thal.tsv"),
    sep = "\t",
    data.table = FALSE
)

###############################################################################



############################# Extract ensembl IDs #############################

message(Sys.time(), " - extracting ensemblIDs")

genes_ensid <- na.omit(unique(sig_genes$ensemblID))
length(genes_ensid)
# [1] 173

jx_ensid <- na.omit(unique(sig_jx$ensemblID))
length(jx_ensid)
# [1] 4592

exon_ensid <- na.omit(unique(sig_exons$ensemblID))
length(exon_ensid)
# [1] 509

###############################################################################



############################# Plot Venn diagrams ##############################

message(Sys.time(), " - plotting Venn diagrams")

plot_venn <- function(list_2plot, category.names, filename, fill) {
    venn.diagram(
        x = list_2plot,
        category.names = category.names,
        filename = filename,
        fill = fill,
        disable.logging = TRUE,
        imagetype = "png",
        units = "in",
        height = 6,
        width = 6,
        total.population = TRUE,
        cat.fontface = "bold",
        cat.default.pos = "outer",
        cex = 1.2,
    )
}

plot_venn(list_2plot = list(genes_ensid, jx_ensid), category.names = c("Gene", "Jx"), filename = here(out_plot, "Venn_gene-jx.png"), fill = c("#f1e7dd", "#eccbd9"))

plot_venn(list_2plot = list(genes_ensid, exon_ensid), category.names = c("Gene", "Exon"), filename = here(out_plot, "Venn_gene-exon.png"), fill = c("#f1e7dd", "#b4e1ff"))

plot_venn(list_2plot = list(exon_ensid, jx_ensid), category.names = c("Exon", "Jx"), filename = here(out_plot, "Venn_exon-jx.png"), fill = c("#b4e1ff", "#eccbd9"))

plot_venn(list_2plot = list(exon_ensid, jx_ensid, genes_ensid), category.names = c("Exon", "Jx", "Gene"), filename = here(out_plot, "Venn_exon-jx-gene.png"), fill = c("#b4e1ff", "#eccbd9", "#f1e7dd"))

###############################################################################
######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
