## June 11, 2023 - Bukola Ajanaku
# Creating official Wallace et al. sce object with cluster info and annotations 
# because they were not in the original sce object and were really hard to find.
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
library("xlsx")

# loading Wallace et al. clean data, updating to version 3 of seurat object so 
# I can just make into a sce object 
wallData <- as.SingleCellExperiment(
  UpdateSeuratObject(
    readRDS(file = here("processed-data", 
                        "09_trans_special_analysis",
                        "Wallace_mouse_data.rds"))
  ))

wallSeu <- readRDS(file = here("processed-data", 
                               "09_trans_special_analysis",
                               "Wallace_mouse_data.rds"))

wall2s <- UpdateSeuratObject(wallSeu)
idents(wallSeu)

############ ADDING SAMPLE INFORMATION #########################################
# realized that the sample data was not in the colData of the Wallace sce object
# last minute
colData(wallData)$Sample <- NA

# grabbing nuclei identifiers
colData(wallData)$Row <- rownames(colData(wallData))

# adding Sample info
colData(wallData)[startsWith(colData(wallData)$Row, "hab_160822" ), ]$Sample <- "Mouse1"
colData(wallData)[startsWith(colData(wallData)$Row, "hab_161102" ), ]$Sample <- "Mouse2"
colData(wallData)[startsWith(colData(wallData)$Row, "hab_161103" ), ]$Sample <- "Mouse3"
colData(wallData)[startsWith(colData(wallData)$Row, "hab_161105" ), ]$Sample <- "Mouse4"

# adding hemispheric information
colData(wallData)$Hemi <- NA
colData(wallData)$Hemi <- ss(x = colData(wallData)$Row, "_", slot = 4)

############ ADDING CLUSTER INFORMATION #########################################
# adding Wallace cell types (found in supplements of eLife paper.)
eLifeCT <- read.csv(file = here("processed-data", "09_trans_special_analysis",
                                "elife_ct.csv"))

# grabbing cluster info
eLife_clusts <- read.xlsx(file = here("processed-data", "09_trans_special_analysis",
                                "elife_clusters.xlsx"), 1, header=TRUE)

names(eLife_clusts)[[1]] <- "gene"

# adding gene names to cell type data
# grabbing Ensembl GRCm38 release 87 information for Wallace data set
ah <- AnnotationHub()
query(ah, "EnsDb.Mmusculus.v87")
org.Mm.eg.db <- ah[["AH53222"]]

# grabbing relevant data from ensdb
mouse_genes <- DataFrame(genes(org.Mm.eg.db))

# making sure wallData has a column for it's symbols
eLifeCT$gene <- rownames(wallData)

# grabbing symbol IDs
eLifeCT$gene_id <- mouse_genes$gene_id[match(eLifeCT$gene, mouse_genes$gene_name)]

table(is.na(eLifeCT$gene_id))
# FALSE  TRUE 
# 22944  2345


# .


