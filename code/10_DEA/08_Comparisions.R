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


######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
