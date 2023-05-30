library("here")
library("SummarizedExperiment")
library("recount")
library("dplyr")
library("ggplot2")
library("cowplot")
library("stringr")
library("RColorBrewer")
library("sessioninfo")



############################# Load rse gene object ############################

load(
    here(
        "preprocessed_data",
        "count_data_bukola",
        "rse_gene_Roche_Habenula_qcAndAnnotated_n69.Rdata"
    ),
    verbose = TRUE
)

lobstr::obj_size(rse_gene)
# 40.63 MB

class(rse_gene)
# [1] "RangedSummarizedExperiment"
# attr(,"package")
# [1] "SummarizedExperiment"

rse_gene
# class: RangedSummarizedExperiment
# dim: 58037 69
# metadata(0):
# assays(1): counts
# rownames(58037): ENSG00000223972.5 ENSG00000227232.5 ...
#   ENSG00000210195.2 ENSG00000210196.2
# rowData names(10): Length gencodeID ... NumTx gencodeTx
# colnames: NULL
# colData names(68): RNum RIN ... Flowcell hasGenotype

dim(colData(rse_gene))
# [1] 69 68

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

unique(colData(rse_gene)$PrimaryDx)
# [1] "Schizo"  "Control"

table(colData(rse_gene)$PrimaryDx)
# Control  Schizo
#      34      35

table(colData(rse_gene)$Sex)
#  M
# 69

table(colData(rse_gene)$Race)
# CAUC
#   69

table(colData(rse_gene)$Flowcell)
# HVYTYBBXX HW252BBXX
#        34        35

###############################################################################



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
