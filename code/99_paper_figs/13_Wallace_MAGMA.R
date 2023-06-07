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
# library("RCurl")
library("AnnotationHub")


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
ah <- AnnotationHub()
query(ah, "EnsDb.Mmusculus.v87")
org.Mm.eg.db <- ah[["AH53222"]]
    # EnsDb for Ensembl:
    #   |Backend: SQLite
    # |Db type: EnsDb
    # |Type of Gene ID: Ensembl Gene ID
    # |Supporting package: ensembldb
    # |Db created by: ensembldb package from Bioconductor
    # |script_version: 0.3.1
    # |Creation time: Fri Jun  9 08:40:26 2017
    # |ensembl_version: 87
    # |ensembl_host: localhost
    # |Organism: mus_musculus
    # |taxonomy_id: 10090
    # |genome_build: GRCm38
    # |DBSCHEMAVERSION: 2.1
    # | No. of genes: 50143.
    # | No. of transcripts: 124168.
    # |Protein data available.

# grabbing relevant data from ensdb
mouse_gene_db <- DataFrame(genes(org.Mm.eg.db))

# making sure wallData has a column for it's symbols
rowData(wallData)$Symbol <- rownames(wallData)

# adding gene_id to rowData of the Wallace sce object
rowData(wallData)$gene_id <- mouse_gene_db$gene_id[match(rowData(wallData)$ID, mouse_gene_db$gene_name)]

table(!is.na(rowData(wallData)$gene_id))
  # FALSE  TRUE 
  # 2345 22944     <- good, 22,944 gene_ids matched

# adding entrez_id to rowData of the Wallace sce object
rowData(wallData)$entrez_id <- as.character(mouse_gene_db$entrezid[match(rowData(wallData)$ID, 
                                                                         mouse_gene_db$symbol)])

table(!is.na(rowData(wallData)$entrez_id))
  # FALSE  TRUE 
  # 3764 21525      <- Not that far off, we're looking good.

######## COMPARING ORGANISMS ###################################################
hom_mm <- hom[hom$Common.Organism.Name == "mouse, laboratory", ]
  # 21829 entries

table(rowData(wallData)$entrez_id %in% hom_mm$EntrezGene.ID)
  # FALSE  TRUE 
  # 6723 18566     <- good, we have 18,566 genes in the shared list

table(rowData(wallData)$ID %in% hom_mm$Symbol)
  # FALSE  TRUE 
  # 8155 17134     <- again, not bad and not that far off!

# adding JAX.geneID by Entrez to wallData
rowData(wallData)$JAX.geneID <- hom_mm$DB.Class.Key[match(rowData(wallData)$entrez_id,
                                                             hom_mm$EntrezGene.ID)]
  # head(rowData(wallData)[ rowData(wallData)$entrez_id != "NULL",])

# Now, we compare the JAX.geneIDs between the two data sets 
length(intersect(rowData(sce)$JAX.geneID,
                 rowData(wallData)$JAX.geneID))  # 15,820
# saving in object
sharedHomologs <- intersect(rowData(sce)$JAX.geneID,
                                      rowData(wallData)$JAX.geneID)
  # [1]       NA 44099149 44108115 44108230 44107312 44102688
  # That first one is NA - rm
sharedHomologs <- sharedHomologs[-1]
  # [1] 44099149 44108115 44108230 44107312 44102688 44096921

# Human not in mouse
length(setdiff(rowData(sce)$JAX.geneID,
               rowData(wallData)$JAX.geneID))  # 1267
# Mouse not in human
length(setdiff(rowData(wallData)$JAX.geneID,
               rowData(sce)$JAX.geneID))  # 2739

# Subset for the shared homologs
sce.mm.sub <- wallData[rowData(wallData)$JAX.geneID %in% sharedHomologs, ]   # 15824
sce.hsap.sub <- sce[rowData(sce)$JAX.geneID %in% sharedHomologs, ]  # 16283
  ## Many are duplicated...

rowData(sce.mm.sub)$ID[duplicated(rowData(sce.mm.sub)$JAX.geneID)]
  # only these duplicates: [1] "Gm37240" "Rsph10b" "Yjefn3"  "Zfp708"  "a"

rowData(sce.hsap.sub)$Symbol[duplicated(rowData(sce.hsap.sub)$JAX.geneID)]
  # total of 464

#### getting rid of the duplicates ############
## Human ===

# first changing the rownames to EnsemblIDs
rownames(sce.hsap.sub) <- rowData(sce.hsap.sub)$ID

duplicatedSet.human <- which(duplicated(rowData(sce.hsap.sub)$JAX.geneID))
genes2compare.human <- list()
gene2keep.human <- character()

for(g in 1:length(duplicatedSet.human)){
  genes2compare.human[[g]] <- rownames(sce.hsap.sub)[rowData(sce.hsap.sub)$JAX.geneID ==
                                        rowData(sce.hsap.sub)$JAX.geneID[duplicatedSet.human[g]]]
  rowmeansmat <- rowMeans(assay(sce.hsap.sub[genes2compare.human[[g]], ], "logcounts"))
  gene2keep.human[g] <- names(rowmeansmat[order(rowmeansmat, decreasing=TRUE)])[1]
}

length(genes2compare.human) # 206
length(unique(gene2keep.human)) # 181


# This is because many 'tested' might have been orthologous,
#   b/tw themselves (i.e. 3+ orthologous genes):
length(unique(rowData(sce.rn.sub)$JAX.geneID[duplicatedSet.rat])) # 181 - good

genesNoCompare.rat <- rownames(sce.rn.sub)[!(rownames(sce.rn.sub) %in% unlist(genes2compare.rat))]

# Finally combine and subset
sce.rn.sub <- sce.rn.sub[c(genesNoCompare.rat, unique(gene2keep.rat)), ]

table(rowData(sce.rn.sub)$JAX.geneID %in% sharedHomologs) # 14007 TRUE
table(duplicated(rowData(sce.rn.sub)$JAX.geneID)) # 14007 FALSE         dope.






# .





