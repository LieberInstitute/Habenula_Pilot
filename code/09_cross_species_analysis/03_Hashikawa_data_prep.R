
library("Seurat")
library("SingleCellExperiment")
library("here")
library("sessioninfo")
library("EnsDb.Hsapiens.v86")
library("org.Hs.eg.db")
library("AnnotationHub")
library("spatialLIBD")
library("Matrix")
library("purrr")

#### load mouse data ####
count_files <- map_chr(c(all = "count.csv", neuron = "count_neuron.csv"), ~here("processed-data", 
                                                                              "09_cross_species_analysis",
                                                                              "Hashikawa_data",
                                                                              .x))

meta_files <- map_chr(c(all = "meta.csv", neuron = "meta_neuron.csv"), ~here("processed-data", 
                                                                     "09_cross_species_analysis",
                                                                     "Hashikawa_data",
                                                                     .x))

all(file.exists(c(count_files, meta_files)))


sce_mouse <- map2(count_files, meta_files, function(count_fn, meta_fn){
  message(Sys.time(), " - Reading: ", count_fn)
  counts <- read.csv(count_fn, row.names = 1)
  
  SingleCellExperiment(
    colData=DataFrame(read.csv(meta_fn, row.names = 1)),
    assays = list(counts = as(counts, "sparseMatrix"))
  )
} )

map(sce_mouse, dim)
# $all
# [1] 17726 11878
# 
# $neuron
# [1] 17726  5558

## empty df
rowData(sce_mouse$all)

map(sce_mouse, ~table(.x$stim, .x$celltype))
# $all
# 
#      Astrocyte1 Astrocyte2 Endothelial Epen Microglia Mural Neuron1 Neuron2 Neuron3 Neuron4 Neuron5 Neuron6 Neuron7
# cntl        969         72          93    6       156    94     429     541     395     461     372     344     181
# stim        637         40          69   36       147   152     631     502     620     326     399     377     130
# 
#      Neuron8 Oligo1 Oligo2 Oligo3 OPC1 OPC2 OPC3
# cntl     225    782    304     77    9  359   73
# stim      61    685    210     50  581  226   57
# 
# $neuron
# 
#      LHb1 LHb2 LHb3 LHb4 LHb5 LHb6 MHb1 MHb2 MHb3 MHb4 MHb5 MHb6
# cntl  329  279  217  214  150  174  315  264  270  148  165  142
# stim  229  178  222  185  210  185  351  398  351  264  229   89


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
# grabbing Ensembl GRCm38 release 87 (?) information for Hashikawa data
ah <- AnnotationHub()
query(ah, "EnsDb.Mmusculus.v87")
org.Mm.eg.db <- ah[["AH53222"]]

# grabbing relevant data from ensdb
mouse_gene_db <- DataFrame(genes(org.Mm.eg.db))
#homolog data
hom_mm <- hom[hom$Common.Organism.Name == "mouse, laboratory", ]


## build row data for sce_mouse
mouse_rd <- DataFrame(Symbol = rownames(sce_mouse$all))
  
# adding gene_id to rowData of the Wallace sce object
mouse_rd$gene_id <- mouse_gene_db$gene_id[match(mouse_rd$Symbol, mouse_gene_db$gene_name)]
table(!is.na(mouse_rd$gene_id))
# FALSE  TRUE 
# 263 17463 
mouse_rd$entrez_id <- as.character(mouse_gene_db$entrezid[match(mouse_rd$Symbol, 
                                                                            mouse_gene_db$symbol)])
table(!is.na(mouse_rd$entrez_id))
# TRUE 
# 17726

#### Compare Organims ####
nrow(hom_mm)
# 21826

table(mouse_rd$entrez_id %in% hom_mm$EntrezGene.ID)
# FALSE  TRUE 
# 2624 15102 

table(mouse_rd$Symbol %in% hom_mm$Symbol)
# FALSE  TRUE 
# 3140 14586 

# adding JAX.geneID by Entrez to sce_mouse
mouse_rd$JAX.geneID <- hom_mm$DB.Class.Key[match(mouse_rd$entrez_id,
                                                          hom_mm$EntrezGene.ID)]

# Find overlapping JAX.geneIDs between the two data sets 
sharedHomologs <- intersect(rowData(sce)$JAX.geneID,
                            mouse_rd$JAX.geneID)
## rm NA
sharedHomologs <- sharedHomologs[!is.na(sharedHomologs)]
length(sharedHomologs) 
# [1] 14473

## add rowData to sce_mouse
rowData(sce_mouse$all) <- mouse_rd
rowData(sce_mouse$neuron) <- mouse_rd

# Subset for the shared homologs
sce_mouse_sub <- map(sce_mouse, ~.x[rowData(.x)$JAX.geneID %in% sharedHomologs, ])
sce_hsap_sub <- sce[rowData(sce)$JAX.geneID %in% sharedHomologs, ]  # 16283

map(sce_mouse_sub, nrow) # [1] 14476
length(sce_hsap_sub) # 14766
## Many are duplicated...

message("Duplicates mm")
length(rowData(sce_mouse_sub$all)$Symbol[duplicated(rowData(sce_mouse_sub$all)$JAX.geneID)])
# 3
rowData(sce_mouse_sub$all)$Symbol[duplicated(rowData(sce_mouse_sub$all)$JAX.geneID)]

message("Duplicates hsap")
length(rowData(sce_hsap_sub)$Symbol[duplicated(rowData(sce_hsap_sub)$JAX.geneID)])
# 293
rowData(sce_hsap_sub)$Symbol[duplicated(rowData(sce_hsap_sub)$JAX.geneID)]

## getting rid of the duplicates
# first changing the rownames to EnsemblIDs
rownames(sce_hsap_sub) <- rowData(sce_hsap_sub)$ID

duplicatedSet.human <- which(duplicated(rowData(sce_hsap_sub)$JAX.geneID))
genes2compare.human <- list()
gene2keep.human <- character()

for(g in 1:length(duplicatedSet.human)){
  genes2compare.human[[g]] <- rownames(sce_hsap_sub)[rowData(sce_hsap_sub)$JAX.geneID ==
                                                       rowData(sce_hsap_sub)$JAX.geneID[duplicatedSet.human[g]]]
  rowmeansmat <- rowMeans(assay(sce_hsap_sub[genes2compare.human[[g]], ], "counts"))
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

sce_mouse_sub <- map(sce_mouse_sub, function(sce_mouse_sub){
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
  
  sce_mouse_sub <- sce_mouse_sub[match(rowData(sce_hsap_sub)$JAX.geneID,
                                       rowData(sce_mouse_sub)$JAX.geneID), ]
  
  # drop count 0 genes
  sce_mouse_sub <- sce_mouse_sub[!rowSums(assay(sce_mouse_sub, "counts"))==0, ]
  
  return(sce_mouse_sub)
})

map_int(sce_mouse_sub, nrow)
#   all neuron 
# 14468  13975 

# subset for corresponding hsap habenula data - match to ALL data
sce_hsap_sub <- sce_hsap_sub[rowData(sce_hsap_sub)$JAX.geneID %in% rowData(sce_mouse_sub$all)$JAX.geneID, ]
length(sce_hsap_sub)
# 14468

all(rowData(sce_hsap_sub)$JAX.geneID == rowData(sce_mouse_sub$all)$JAX.geneID) ## TRUE

## save sce_mouse_sub
save(sce_mouse_sub, file = here("processed-data", 
                                "09_cross_species_analysis",
                                "Hashikawa_data",
                                "sce_mouse_habenula.Rdata"))

## pseudobulk + register data

mouse_modeling_results <- map2(sce_mouse_sub, names(sce_mouse_sub), 
                              ~registration_wrapper( .x,
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
                                                                                paste0("sce_mouse_habenula_pb-",.y,".rds"))
                              ))


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

##subset hsap to just Habenula neurons

sce_hsap_sub <- sce_hsap_sub[,grep("Hb", sce_hsap_sub$final_Annotations)]
ncol(sce_hsap_sub)

hsap_modeling_results_neuron <- registration_wrapper(
  sce_hsap_sub,
  var_registration = "final_Annotations",
  var_sample_id = "Sample",
  covars = NULL,
  gene_ensembl = "JAX.geneID",
  gene_name = "Symbol",
  suffix = "",
  min_ncells = 10)


hsap_modeling_results = list(all = hsap_modeling_results, neuron = hsap_modeling_results_neuron)

save(hsap_modeling_results, mouse_modeling_results, file = here("processed-data", 
                                                                    "09_cross_species_analysis",
                                                                    "Hashikawa_homolog_modeling_results.Rdata"))

# sgejobs::job_single('03_Hashikawa_data_prep', create_shell = TRUE, memory = '30G', command = "Rscript 03_Hashikawa_data_prep.R")


## Reproducibility information
print("Reproducibility information:")
options(width = 120)
session_info()

