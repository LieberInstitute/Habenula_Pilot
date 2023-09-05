
library("Seurat")
library("SingleCellExperiment")
library("here")
library("sessioninfo")
library("EnsDb.Hsapiens.v86")
library("org.Hs.eg.db")
library("AnnotationHub")
library("spatialLIBD")
library(Matrix)

# 
# sce_mouse <- as.SingleCellExperiment(
#   UpdateSeuratObject(
#     readRDS(file = here("processed-data", 
#                         "09_cross_species_analysis",
#                         "Hashikawa_data",
#                         "Habenula_Seurat_object.rds"))
#   ))
# 
# dim(sce_mouse)
# # dim: 2000 11878

mouse_pd <- read.csv(here("processed-data", 
                          "09_cross_species_analysis",
                          "Hashikawa_data",
                          "meta.csv"
                          ), row.names = 1)

head(mouse_pd)

counts <- Matrix(read.csv(here("processed-data", 
                                "09_cross_species_analysis",
                                "Hashikawa_data",
                                "count.csv"), row.names = 1), 
                 sparse = TRUE)

ncol(counts) == nrow(mouse_pd)
all(colnames(counts) == rownames(mouse_pd))

sce_mouse <- SingleCellExperiment(colData=DataFrame(mouse_pd),
                                  assays = list(counts = counts))
# dim: 17726 11878
rm(counts)

## empty df
rowData(sce_mouse)

table(sce_mouse$stim, sce_mouse$celltype)
#      Astrocyte1 Astrocyte2 Endothelial Epen Microglia Mural Neuron1 Neuron2 Neuron3 Neuron4 Neuron5 Neuron6 Neuron7
# cntl        969         72          93    6       156    94     429     541     395     461     372     344     181
# stim        637         40          69   36       147   152     631     502     620     326     399     377     130
# 
#      Neuron8 Oligo1 Oligo2 Oligo3 OPC1 OPC2 OPC3
# cntl     225    782    304     77    9  359   73
# stim      61    685    210     50  581  226   57


#### Human Data ####
# loading our final sce object
load(file = here("processed-data","sce_objects", "sce_Habenula_Pilot.Rdata"), verbose = TRUE)

hs.entrezIds <- mapIds(org.Hs.eg.db, keys = rowData(sce)$ID, 
                       column = "ENTREZID", keytype="ENSEMBL")
table(!is.na(hs.entrezIds))
# FALSE  TRUE 
# 11156 22692 

# storing genes without entrez IDs 
withoutEntrez <- names(hs.entrezIds)[is.na(hs.entrezIds)]
# saving the genes without the EntezIDs elsewhere with their Symbol identities
table(rowData(sce)[rowData(sce)$ID %in% withoutEntrez, ]$ID == withoutEntrez)
names(withoutEntrez) <- rowData(sce)[rowData(sce)$ID %in% withoutEntrez, ]$Symbol

# Add to rowData
rowData(sce) <- cbind(rowData(sce), hs.entrezIds)

# JAX annotation info
hom = read.delim("http://www.informatics.jax.org/downloads/reports/HOM_AllOrganism.rpt",
                 as.is=TRUE)
hom_hs <- hom[hom$Common.Organism.Name == "human", ]
dim(hom_hs)

# adding JAX annotations to our sce metadata by the entrez ID
rowData(sce)$JAX.geneID <- hom_hs$DB.Class.Key[match(rowData(sce)$hs.entrezIds,
                                                     hom_hs$EntrezGene.ID)]

#### Mouse Data ####
# grabbing Ensembl GRCm38 release 87 information for Wallace data set
ah <- AnnotationHub()
query(ah, "EnsDb.Mmusculus.v87")
org.Mm.eg.db <- ah[["AH53222"]]

# grabbing relevant data from ensdb
mouse_gene_db <- DataFrame(genes(org.Mm.eg.db))

# making sure sce_mouse has a column for it's symbols
rowData(sce_mouse)$Symbol <- rownames(sce_mouse)

# adding gene_id to rowData of the Wallace sce object
rowData(sce_mouse)$gene_id <- mouse_gene_db$gene_id[match(rowData(sce_mouse)$Symbol, mouse_gene_db$gene_name)]

table(!is.na(rowData(sce_mouse)$gene_id))
# FALSE  TRUE 
# 263 17463

hom_mm <- hom[hom$Common.Organism.Name == "mouse, laboratory", ]

rowData(sce_mouse)$entrez_id <- as.character(mouse_gene_db$entrezid[match(rowData(sce_mouse)$Symbol, 
                                                                         mouse_gene_db$symbol)])
table(!is.na(rowData(sce_mouse)$entrez_id))

#### Compare Organims ####
nrow(hom_mm)
# 21826

table(rowData(sce_mouse)$entrez_id %in% hom_mm$EntrezGene.ID)
# FALSE  TRUE 
# 2624 15102 

table(rowData(sce_mouse)$Symbol %in% hom_mm$Symbol)
# FALSE  TRUE 
# 3140 14586 

# adding JAX.geneID by Entrez to sce_mouse
rowData(sce_mouse)$JAX.geneID <- hom_mm$DB.Class.Key[match(rowData(sce_mouse)$entrez_id,
                                                          hom_mm$EntrezGene.ID)]

# Now, we compare the JAX.geneIDs between the two data sets 
length(intersect(rowData(sce)$JAX.geneID,
                 rowData(sce_mouse)$JAX.geneID))
# [1] 14474

# saving in object
sharedHomologs <- intersect(rowData(sce)$JAX.geneID,
                            rowData(sce_mouse)$JAX.geneID)
sharedHomologs <- sharedHomologs[!is.na(sharedHomologs)]
length(sharedHomologs) 
# [1] 14473


# Subset for the shared homologs
sce_mouse_sub <- sce_mouse[rowData(sce_mouse)$JAX.geneID %in% sharedHomologs, ]
sce_hsap_sub <- sce[rowData(sce)$JAX.geneID %in% sharedHomologs, ]  # 16283

length(sce_mouse_sub) # [1] 14476
length(sce_hsap_sub) # 14766
## Many are duplicated...

length(rowData(sce_mouse_sub)$Symbol[duplicated(rowData(sce_mouse_sub)$JAX.geneID)])
# 3

length(rowData(sce_hsap_sub)$Symbol[duplicated(rowData(sce_hsap_sub)$JAX.geneID)])
# 293

## getting rid of the duplicates
# first changing the rownames to EnsemblIDs
rownames(sce_hsap_sub) <- rowData(sce_hsap_sub)$ID

duplicatedSet.human <- which(duplicated(rowData(sce_hsap_sub)$JAX.geneID))
genes2compare.human <- list()
gene2keep.human <- character()

for(g in 1:length(duplicatedSet.human)){
  genes2compare.human[[g]] <- rownames(sce_hsap_sub)[rowData(sce_hsap_sub)$JAX.geneID ==
                                                       rowData(sce_hsap_sub)$JAX.geneID[duplicatedSet.human[g]]]
  rowmeansmat <- rowMeans(assay(sce_hsap_sub[genes2compare.human[[g]], ], "logcounts"))
  gene2keep.human[g] <- names(rowmeansmat[order(rowmeansmat, decreasing=TRUE)])[1]
}

length(genes2compare.human)
# [1] 293
length(unique(gene2keep.human))
# [1] 216

# This is because many 'tested' might have been orthologous,
#   b/tw themselves (i.e. 3+ orthologous genes):
length(unique(rowData(sce_hsap_sub)$JAX.geneID[duplicatedSet.human])) # [1] 216

genesNoCompare.human <- rownames(sce_hsap_sub)[!(rownames(sce_hsap_sub) 
                                                 %in% unlist(genes2compare.human))]

# Finally combine and subset
sce_hsap_sub <- sce_hsap_sub[c(genesNoCompare.human, unique(gene2keep.human)), ]

table(rowData(sce_hsap_sub)$JAX.geneID %in% sharedHomologs)
# TRUE 
# 14473
any(duplicated(rowData(sce_hsap_sub)$JAX.geneID))
# [1] FALSE

## Mouse duplicates
# make sure rownames are EnsemblIDs
rownames(sce_mouse_sub) <- rowData(sce_mouse_sub)$gene_id

duplicatedSet.mouse <- which(duplicated(rowData(sce_mouse_sub)$JAX.geneID))
genes2compare.mouse <- list()
gene2keep.mouse <- character()

for(g in 1:length(duplicatedSet.mouse)){
  genes2compare.mouse[[g]] <- rownames(sce_mouse_sub)[rowData(sce_mouse_sub)$JAX.geneID ==
                                                     rowData(sce_mouse_sub)$JAX.geneID[duplicatedSet.mouse[g]]]
  rowmeansmat <- rowMeans(assay(sce_mouse_sub[genes2compare.mouse[[g]], ], "counts"))
  gene2keep.mouse[g] <- names(rowmeansmat[order(rowmeansmat, decreasing=TRUE)])[1]
}

length(unique(rowData(sce_mouse_sub)$JAX.geneID[duplicatedSet.mouse])) # 3

genesNoCompare.mouse <- rownames(sce_mouse_sub)[!(rownames(sce_mouse_sub) %in% unlist(genes2compare.mouse))]

# Finally combine and subset
sce_mouse_sub <- sce_mouse_sub[c(genesNoCompare.mouse, unique(gene2keep.mouse)), ]

#### Save Mouse Data ####
sce_mouse_sub <- sce_mouse_sub[match(rowData(sce_hsap_sub)$JAX.geneID,
                               rowData(sce_mouse_sub)$JAX.geneID), ]

dim(sce_mouse_sub)
# [1] 14473 11878

# drop count 0 genes
sce_mouse_sub <- sce_mouse_sub[!rowSums(assay(sce_mouse_sub, "counts"))==0, ]
dim(sce_mouse_sub)
# [1] 14468 11878

# subset for corresponding hsap habenula data
sce_hsap_sub <- sce_hsap_sub[rowData(sce_hsap_sub)$JAX.geneID %in% rowData(sce_mouse_sub)$JAX.geneID, ]
length(sce_hsap_sub)
# 14468

all(rowData(sce_hsap_sub)$JAX.geneID == rowData(sce_mouse_sub)$JAX.geneID) ## TRUE

## save sce_mouse_sub
save(sce_mouse_sub, file = here("processed-data", 
                                "09_cross_species_analysis",
                                "Hashikawa_data",
                                "sce_mouse_habenula.Rdata"))

## pseudobulk + register data

colData(sce_mouse_sub)

counts(sce_mouse_sub)

mouse_modeling_results <- registration_wrapper(
  sce_mouse_sub,
  var_registration = "celltype",
  var_sample_id = "stim",
  covars = NULL,
  gene_ensembl = "JAX.geneID",
  gene_name = "Symbol",
  suffix = "",
  min_ncells = 10,
  pseudobulk_rds_file = here("processed-data", 
                             "09_cross_species_analysis",
                             "Hashikawa_data",
                             "sce_mouse_habenula_pb.rds")
)

mouse_modeling_results 

hsap_modeling_results <- registration_wrapper(
  sce_hsap_sub,
  var_registration = "final_Annotations",
  var_sample_id = "Sample",
  covars = NULL,
  gene_ensembl = "JAX.geneID",
  gene_name = "Symbol",
  suffix = "",
  min_ncells = 10,
  pseudobulk_rds_file = here("processed-data", 
                             "09_cross_species_analysis",
                             "Hashikawa_data",
                             "sce_habenula_pb.rds")
)


save(hsap_modeling_results, mouse_modeling_results, file = here("processed-data", 
                                                                    "09_cross_species_analysis",
                                                                    "modeling_results.Rdata"))

## human habenula t-stats
registration_t_stats <- hsap_modeling_results$enrichment[, grep("^t_stat", colnames(hsap_modeling_results$enrichment))]
colnames(registration_t_stats) <- gsub("^t_stat_", "", colnames(registration_t_stats))
rownames(registration_t_stats) <- hsap_modeling_results$enrichment$ensembl

cor_species <- layer_stat_cor(
  modeling_results = mouse_modeling_results,
  stats = registration_t_stats,
  model_type = "enrichment",
  top_n = 100
)

pdf(here("plots", "09_cross_species_analysis","cross_species.pdf"))
layer_stat_cor_plot(cor_species, max = max(cor_species))
dev.off()

## Reproducibility information
print("Reproducibility information:")
options(width = 120)
session_info()

