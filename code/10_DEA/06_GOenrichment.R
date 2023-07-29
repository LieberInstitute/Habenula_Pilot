library("here")
library("data.table")
library("dplyr")
library("clusterProfiler")
library("org.Hs.eg.db")
library("ggplot2")
library("sessioninfo")

out_plot <- here("plots", "10_DEA", "05_GOenrichment")
out_data <- here("processed-data", "10_DEA", "05_GOenrichment")



####################### Load tsv with all genes from DEA ######################

DE_all_files <- list.files(
    here(
        "processed-data",
        "10_DEA",
        "04_DEA"
    ),
    pattern = "DEA_All*",
    full.names = TRUE
)
DE_all_files <- DE_all_files[1:3]

DE_all <- lapply(DE_all_files, fread, data.table = FALSE)

DE_sig_files <- list.files(
    here(
        "processed-data",
        "10_DEA",
        "04_DEA"
    ),
    pattern = "DEA_Sig*",
    full.names = TRUE
)
DE_sig_files <- DE_sig_files[1:3]

DE_sig <- lapply(DE_sig_files, fread, data.table = FALSE)

###############################################################################



######################### Prepare data fo GO analysis #########################

## Converting ensembl IDs to Entrez IDs in case some of my data.frames does not include that column
DE_all <- lapply(DE_all, function(x) {
    entrez <- mapIds(org.Hs.eg.db, keys = x$ensemblID, keytype = "ENSEMBL", column = "ENTREZID")
    entrez[is.na(names(entrez))] <- NA
    entrez <- unlist(entrez)
    x$EntrezID <- entrez
    return(x)
})

DE_sig <- lapply(DE_sig, function(x) {
    entrez <- mapIds(org.Hs.eg.db, keys = x$ensemblID, keytype = "ENSEMBL", column = "ENTREZID")
    entrez[is.na(names(entrez))] <- NA
    entrez <- unlist(entrez)
    x$EntrezID <- entrez
    return(x)
})

## Split in: all, down, up and using just unique EntrezIDs
sigFeat <- lapply(DE_sig, function(x) {
    sigGene <- c(list(x$EntrezID), split(x$EntrezID, sign(x$logFC)))
    sigGene <- lapply(sigGene, function(x) x[!is.na(x)])
    sigGene <- lapply(sigGene, unique)
    names(sigGene) <- c("all", "down", "up")
    return(sigGene)
})

## Using just unique EntrezIDs
allFeat <- lapply(DE_all, function(x) {
    allGene <- as.character(x$EntrezID)
    allGene <- allGene[!is.na(allGene)]
    allGene <- unique(allGene)
    return(allGene)
})

###############################################################################



##################### Run GO and KEGG enrichment analysis #####################

go <- compareCluster(sigGene,
    fun = "enrichGO",
    universe = geneUniverse,
    OrgDb = org.Hs.eg.db,
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.15,
    #qvalueCutoff = 0.05,
    readable = TRUE
)

kegg <- compareCluster(sigGene,
    fun = "enrichKEGG",
    universe = geneUniverse,
    organism = "human",
    pvalueCutoff = 0.15,
    qvalueCutoff = 0.05
)

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
