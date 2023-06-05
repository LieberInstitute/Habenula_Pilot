## June 5, 2023 - Bukola Ajanaku
# Working on MAGMA trans-special analysis of our sce object against the 
# Wallace et al. 2019 paper mouse data.
# qrsh -l mem_free=50G,h_vmem=50G

library("SingleCellExperiment")
library("here")
library("sessioninfo")
library("jaffelab")
library("Seurat")

# loading our final sce object
load(file = here("processed-data", "99_paper_figs", "sce_objects", 
                 "official_final_sce.RDATA"), verbose = TRUE)
  # sce 

# checking cluster data 
table(sce$final_Annotations)

# loading Wallace et al. clean data
wallData <- readRDS(file = here("processed-data", "99_paper_figs", "MAGMA",
                 "Wallace_mouse_data.rds"))

test <- as.data.frame(wallData)


# sourcing official color palette 
source(file = here("code", "99_paper_figs", "source_colors.R"))
  # bulk_colors and sn_colors

# adaptation of Matt's code:
# Add EntrezID for human
hs.entrezIds <- mapIds(org.Hs.eg.db, keys=rowData(sce.nac)$gene_id, 
                       column="ENTREZID", keytype="ENSEMBL")
# "'select()' returned 1:many mapping between keys and columns"
table(!is.na(hs.entrezIds))
# 21,191 valid entries (remember this is already subsetted for those non-zero genes only)
withoutEntrez <- names(hs.entrezIds)[is.na(hs.entrezIds)]
# Store those somewhere, maybe for later reference
table(rowData(sce.nac)[rowData(sce.nac)$gene_id %in% withoutEntrez, ]$gene_id == withoutEntrez)
names(withoutEntrez) <- rowData(sce.nac)[rowData(sce.nac)$gene_id %in% withoutEntrez, ]$gene_name


# Add to rowData
rowData(sce.nac) <- cbind(rowData(sce.nac), hs.entrezIds)
