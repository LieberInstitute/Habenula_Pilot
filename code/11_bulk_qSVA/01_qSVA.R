library("here")
library("SummarizedExperiment")
library("qsvaR")
library("sessioninfo")



############################### Load rse objects ##############################

load(
    here(
        "processed-data",
        "02_bulk_qc",
        "count_data_bukola",
        "rse_gene_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata"
    ),
    verbose = TRUE
)
# Loading objects:
#   rse_gene

load(
    here(
        "processed-data",
        "02_bulk_qc",
        "count_data_bukola",
        "rse_tx_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata"
    ),
    verbose = TRUE
)
# Loading objects:
#   rse_tx

rse_gene
# class: RangedSummarizedExperiment
# dim: 22756 69
# metadata(0):
# assays(2): counts logcounts
# rownames(22756): ENSG00000227232.5 ENSG00000278267.1 ...
#   ENSG00000210195.2 ENSG00000210196.2
# rowData names(11): Length gencodeID ... gencodeTx MGI_Symbol
# colnames(69): R18346 R18347 ... R18423 R18424
# colData names(76): RNum RIN ... subsets_Ribo_detected
#   subsets_Ribo_percent

dim(rse_gene)
# [1] 22756    69

names(colData(rse_gene))
#  [1] "RNum"                           "RIN"
#  [3] "Brain.Region"                   "BrNum"
#  [5] "AgeDeath"                       "Sex"
#  [7] "Race"                           "PrimaryDx"
#  [9] "FQCbasicStats"                  "perBaseQual"
# [11] "perTileQual"                    "perSeqQual"
# [13] "perBaseContent"                 "GCcontent"
# [15] "Ncontent"                       "SeqLengthDist"
# [17] "SeqDuplication"                 "OverrepSeqs"
# [19] "AdapterContent"                 "KmerContent"
# [21] "SeqLength_R1"                   "percentGC_R1"
# [23] "phred20.21_R1"                  "phred48.49_R1"
# [25] "phred76.77_R1"                  "phred100.101_R1"
# [27] "phredGT30_R1"                   "phredGT35_R1"
# [29] "Adapter50.51_R1"                "Adapter70.71_R1"
# [31] "Adapter88.89_R1"                "SeqLength_R2"
# [33] "percentGC_R2"                   "phred20.21_R2"
# [35] "phred48.49_R2"                  "phred76.77_R2"
# [37] "phred100.101_R2"                "phredGT30_R2"
# [39] "phredGT35_R2"                   "Adapter50.51_R2"
# [41] "Adapter70.71_R2"                "Adapter88.89_R2"
# [43] "ERCCsumLogErr"                  "bamFile"
# [45] "trimmed"                        "numReads"
# [47] "numMapped"                      "numUnmapped"
# [49] "overallMapRate"                 "concordMapRate"
# [51] "totalMapped"                    "mitoMapped"
# [53] "mitoRate"                       "totalAssignedGene"
# [55] "gene_Assigned"                  "gene_Unassigned_Ambiguity"
# [57] "gene_Unassigned_MultiMapping"   "gene_Unassigned_NoFeatures"
# [59] "gene_Unassigned_Unmapped"       "gene_Unassigned_MappingQuality"
# [61] "gene_Unassigned_FragmentLength" "gene_Unassigned_Chimera"
# [63] "gene_Unassigned_Secondary"      "gene_Unassigned_Nonjunction"
# [65] "gene_Unassigned_Duplicate"      "rRNA_rate"
# [67] "Flowcell"                       "hasGenotype"
# [69] "sum"                            "detected"
# [71] "subsets_Mito_sum"               "subsets_Mito_detected"
# [73] "subsets_Mito_percent"           "subsets_Ribo_sum"
# [75] "subsets_Ribo_detected"          "subsets_Ribo_percent"


rse_tx
# class: RangedSummarizedExperiment
# dim: 82434 69
# metadata(0):
# assays(2): tpm logcounts
# rownames(82434): ENST00000488147.1 ENST00000461467.1 ...
#   ENST00000387460.2 ENST00000387461.2
# rowData names(22): source type ... protein_id ccdsid
# colnames: NULL
# colData names(68): RNum RIN ... Flowcell hasGenotype

dim(rse_tx)
# [1] 82434    69

names(colData(rse_tx))
#  [1] "RNum"                           "RIN"
#  [3] "Brain.Region"                   "BrNum"
#  [5] "AgeDeath"                       "Sex"
#  [7] "Race"                           "PrimaryDx"
#  [9] "FQCbasicStats"                  "perBaseQual"
# [11] "perTileQual"                    "perSeqQual"
# [13] "perBaseContent"                 "GCcontent"
# [15] "Ncontent"                       "SeqLengthDist"
# [17] "SeqDuplication"                 "OverrepSeqs"
# [19] "AdapterContent"                 "KmerContent"
# [21] "SeqLength_R1"                   "percentGC_R1"
# [23] "phred20.21_R1"                  "phred48.49_R1"
# [25] "phred76.77_R1"                  "phred100.101_R1"
# [27] "phredGT30_R1"                   "phredGT35_R1"
# [29] "Adapter50.51_R1"                "Adapter70.71_R1"
# [31] "Adapter88.89_R1"                "SeqLength_R2"
# [33] "percentGC_R2"                   "phred20.21_R2"
# [35] "phred48.49_R2"                  "phred76.77_R2"
# [37] "phred100.101_R2"                "phredGT30_R2"
# [39] "phredGT35_R2"                   "Adapter50.51_R2"
# [41] "Adapter70.71_R2"                "Adapter88.89_R2"
# [43] "ERCCsumLogErr"                  "bamFile"
# [45] "trimmed"                        "numReads"
# [47] "numMapped"                      "numUnmapped"
# [49] "overallMapRate"                 "concordMapRate"
# [51] "totalMapped"                    "mitoMapped"
# [53] "mitoRate"                       "totalAssignedGene"
# [55] "gene_Assigned"                  "gene_Unassigned_Ambiguity"
# [57] "gene_Unassigned_MultiMapping"   "gene_Unassigned_NoFeatures"
# [59] "gene_Unassigned_Unmapped"       "gene_Unassigned_MappingQuality"
# [61] "gene_Unassigned_FragmentLength" "gene_Unassigned_Chimera"
# [63] "gene_Unassigned_Secondary"      "gene_Unassigned_Nonjunction"
# [65] "gene_Unassigned_Duplicate"      "rRNA_rate"
# [67] "Flowcell"                       "hasGenotype"

###############################################################################



##################### Load deconvolution and SNP PCs data #####################

## Load SNP PCs
snpPCs <- read.table(
    here(
        "processed-data",
        "08_bulk_snpPC",
        "habenula_genotypes_filt.snpPCs.tab"
    ),
    header = TRUE
)

## Load deconvolution results
load(
    here(
        "processed-data",
        "06_deconvolution",
        "run_Bisque",
        "est_prop_split_Hb_annotations.RDATA"
    ),
    verbose = TRUE
)

###############################################################################



####################### Prepare data to add to the model ######################

## Add some metrics to colData
colData(rse_gene)$library_size <- apply(assay(rse_gene), 2, sum)
colData(rse_gene)$log10_library_size <- log10(colData(rse_gene)$library_size)
colData(rse_gene)$detected_num_genes <- apply(assay(rse_gene), 2, function(x) {
    length(x[which(x > 0)])
})
colData(rse_gene)$abs_ERCCsumLogErr <- abs(colData(rse_gene)$ERCCsumLogErr)

## Add total proportion of Hb (LHb + MHb) and total Thal (Excit.Thal + Inhib.Thal)
est_prop <- t(est_prop$bulk.props)

est_prop <- cbind(est_prop, est_prop[, colnames(est_prop) == "LHb"] + est_prop[, colnames(est_prop) == "MHb"])
est_prop <- cbind(est_prop, est_prop[, colnames(est_prop) == "Excit.Thal"] + est_prop[, colnames(est_prop) == "Inhib.Thal"])

colnames(est_prop)[10:11] <- c("tot.Hb", "tot.Thal")

## Delete sample Br5572/R18424
rse_tx <- rse_tx[, rse_tx$BrNum != "Br5572"]
rse_gene <- rse_gene[, rse_gene$BrNum != "Br5572"]
est_prop <- est_prop[rownames(est_prop) != "R18424", ]

###############################################################################



############## Add cell proportions and SNP PCs data to colData ###############

## Add cell proportions
rownames(est_prop) == rownames(colData(rse_gene)) ## I'm checking if samples are in the same order
colData(rse_gene) <- cbind(colData(rse_gene), est_prop)

## Add SNP PCs
# merged_col <- merge(colData(rse_gene), as.data.frame(snpPCs), by = "BrNum", sort = FALSE)
# stopifnot(identical(rse_gene$RNum, merged_col$RNum))
# colData(rse_gene) <- merged_col
# colnames(rse_gene) <- colData(rse_gene)$RNum

## Copy colData() from rse_gene to rse_tx
colData(rse_tx) <- colData(rse_gene)

###############################################################################



################################## Set model ##################################

mod <- model.matrix(~ PrimaryDx + AgeDeath + Flowcell + mitoRate + rRNA_rate + totalAssignedGene + RIN + abs_ERCCsumLogErr + tot.Hb + tot.Thal,
    data = colData(rse_tx)
)

###############################################################################



################################### Run qSVA ##################################

set.seed(20230703)
qsva_pcs_standard <- qsvaR::qSVA(rse_tx, type = "standard", mod = mod, assayname = "tpm")
dim(qsva_pcs_standard)
# [1] 68  8

set.seed(20230703)
qsva_pcs_cc <- qSVA(rse_tx, type = "cell_component", mod = mod, assayname = "tpm")
dim(qsva_pcs_cc)
# [1] 68  10

###############################################################################



##################### Explore differences between types #######################

rse_cellcomp <- getDegTx(rse_tx, type = "cell_component")
dim(rse_cellcomp)
# [1] 2938   68

rse_stand <- getDegTx(rse_tx, type = "standard")
dim(rse_stand)
# [1] 1772   68

## Aparently the standard just has less genes
length(intersect(rownames(rse_stand), rownames(rse_cellcomp)))
# [1] 1772
length(union(rownames(rse_stand), rownames(rse_cellcomp)))
# [1] 2938

###############################################################################



########################## Write table with qsva data #########################

write.table(qsva_pcs_standard,
    here(
        "processed-data",
        "11_bulk_qSVA",
        "qSVA.tsv"
    ),
    sep = "\t",
    quote = FALSE
)

###############################################################################



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
