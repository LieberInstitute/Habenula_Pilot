
library("SummarizedExperiment")
library("here")
library("sessioninfo")

## dirs
data_dir <- here("processed-data", "03_bulk_pca", "02_multiregion_PCA")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

# plot_dir <- here("plots", "03_bulk_pca", "02_multiregion_PCA")
# if(!dir.exists(plot_dir)) dir.create(plot, recursive = TRUE)

#### load data ####
## retreiv colData final rse 
load(here( "processed-data","rse_objects","rse_gene_Habenula_Pilot.rda"),verbose = TRUE)
pd <- colData(rse_gene)

## unfiltered rse_gene
load(here("preprocessed_data", "rse_gene_Roche_Habenula_PairedEnd_n73.Rdata")) # gene info
rse_gene_Hb <- rse_gene
dim(rse_gene_Hb)
# [1] 58037    73

rse_gene_Hb <- rse_gene_Hb[,rownames(pd)]
colData(rse_gene_Hb) <- pd

## control only
rse_gene_Hb <- rse_gene_Hb[,rse_gene_Hb$PrimaryDx == "Control"]
dim(rse_gene_Hb)
# [1] 58037    33

## other region data 
load("/dcs04/lieber/lcolladotor/dbDev_LIBD001/subsets/for_Louise_230628_nonHabenula/rse_gene_ribozero_nonHabenula_n784.rda", verbose = TRUE)
dim(rse_gene)
# [1] 58037   784

table(rse_gene$Region)
# Amygdala      BLA       CA  Caudate     dACC       DG    DLPFC    HIPPO      MeA     mPFC     sACC 
#      140       54       11       82       55       41      121       26       55       57      142 

table(rse_gene$Dataset)
# Astellas_DG   BrainSeq_Phase2_DLPFC   BrainSeq_Phase2_HIPPO BrainSeq_Phase3_Caudate     BrainSeq_Phase4and5 
#          30                      66                      12                      82                      32 
# psychENCODE_Mood         PTSD_BrainOmics                 VA_PTSD 
#              268                      75                     219 


colnames(colData(rse_gene))
# [1] "SAMPLE_ID"         "RNum"              "RIN"               "Region"            "Dataset"          
# [6] "BrNum"             "Dx"                "Age"               "Sex"               "Race"             
# [11] "Protocol"          "numReads"          "numMapped"         "numUnmapped"       "mitoMapped"       
# [16] "totalMapped"       "overallMapRate"    "concordMapRate"    "mitoRate"          "rRNA_rate"        
# [21] "totalAssignedGene" "bamFile" 

colnames(colData(rse_gene_Hb))

## Add missing cols to Hb data

colnames(colData(rse_gene_Hb))[colnames(colData(rse_gene_Hb)) %in% colnames(colData(rse_gene))]

colnames(colData(rse_gene))[!colnames(colData(rse_gene)) %in% colnames(colData(rse_gene_Hb))]
# [1] "SAMPLE_ID" "Region"    "Dataset"   "Dx"        "Age"       "Protocol" 


