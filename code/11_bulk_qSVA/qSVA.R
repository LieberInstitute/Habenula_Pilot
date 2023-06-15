library("here")
library("qsvaR")
library("sessioninfo")



############################ Load rse gene objects ############################

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


######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
