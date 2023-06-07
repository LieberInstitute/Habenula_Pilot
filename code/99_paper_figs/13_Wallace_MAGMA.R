## June 5, 2023 - Bukola Ajanaku
# Working on MAGMA trans-special analysis of our sce object against the 
# Wallace et al. 2019 paper mouse data.
# qrsh -l mem_free=50G,h_vmem=50G

library("SingleCellExperiment")
library("here")
library("sessioninfo")
library("jaffelab")
library("Seurat")
library("scater")
library("EnsDb.Hsapiens.v86")
#BiocManager::install("org.Rn.eg.db")
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("Mus.musculus")
# library("wget")
# library("BSgenome.Mmusculus.UCSC.mm10")
library("RCurl")

# loading our final sce object
load(file = here("processed-data", "99_paper_figs", "sce_objects", 
                 "official_final_sce.RDATA"), verbose = TRUE)
  # sce 

# checking cluster data 
table(sce$final_Annotations)

# loading Wallace et al. clean data, updating to version 3 of seurat object so 
# I can just make into a sce object 
wallData <- as.SingleCellExperiment(
              UpdateSeuratObject(
                readRDS(file = here("processed-data", "99_paper_figs", "MAGMA",
                "Wallace_mouse_data.rds"))
              ))
  
# sourcing official color palette 
source(file = here("code", "99_paper_figs", "source_colors.R"))
  # bulk_colors and sn_colors

#### Adaptation of Matt's code #################################################
####### HUMAN ##################################################################
# Add EntrezID for human genes in our final sce 
hs.entrezIds <- mapIds(org.Hs.eg.db, keys = rowData(sce)$ID, 
                       column = "ENTREZID", keytype="ENSEMBL")
  # "'select()' returned 1:many mapping between keys and columns"

table(!is.na(hs.entrezIds))
  # FALSE  TRUE 
  # 11126 22722

# adding infor to metaData of our sce object 

# storing genes without entrez IDs 
withoutEntrez <- names(hs.entrezIds)[is.na(hs.entrezIds)]
# saving the genes without the EntezIDs elsewhere with their Symbol identities
table(rowData(sce)[rowData(sce)$ID %in% withoutEntrez, ]$ID == withoutEntrez)
names(withoutEntrez) <- rowData(sce)[rowData(sce)$ID %in% withoutEntrez, ]$Symbol

# Add to rowData
rowData(sce) <- cbind(rowData(sce), hs.entrezIds)

# Bring in 'DB.Class.Key' for human ===
# JAX annotation info
hom = read.delim("http://www.informatics.jax.org/downloads/reports/HOM_AllOrganism.rpt",
                 as.is=TRUE)
hom_hs <- hom[hom$Common.Organism.Name == "human", ]
dim(hom_hs)
# [1] 24609    12      <- 24,609 entries

table(rowData(sce)$hs.entrezIds %in% hom_hs$EntrezGene.ID)
# 17,632 trues
table(rowData(sce)$Symbol %in% hom_hs$Symbol)
# 17,249 trues - very minor difference which is good

# adding JAX annotations to our sce metadata by the entrez ID
rowData(sce)$JAX.geneID <- hom_hs$DB.Class.Key[match(rowData(sce)$hs.entrezIds,
                                                     hom_hs$EntrezGene.ID)]

####### MOUSE ##################################################################
# grabbing Ensembl GRCm38 release 87 information for Wallace data set


rowData(wallData)$ID <- rownames(wallData)

# grab the Symbol information 
Mm.Symbol <- mapIds(org.Mm.eg.db, keys=rowData(wallData)$ID, 
                    column="ENSEMBL", keytype="SYMBOL")

table(!is.na(Mm.Symbol))
# FALSE  TRUE 
# 6515 18774 

# Add Symbols to rowData
rowData(wallData) <- cbind(rowData(wallData), Mm.Symbol)

# finding entrez IDs 
Mm.entrezIds <- mapIds(BSgenome.Mmusculus.UCSC.mm10, keys=rowData(wallData)$ID, 
                     column="ENTREZID", keytype="SYMBOL")

table(!is.na(Mm.entrezIds))










# 





