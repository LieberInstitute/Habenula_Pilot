library(SingleCellExperiment)
library(scRNAseq)
library(batchelor)
library(EnsDb.Hsapiens.v86)
library(scater)
library(scran)
library(scry)
library(uwot)
library(DropletUtils)
library(jaffelab)
library(Rtsne)
library(gridExtra)


### Palette taken from `scater`
tableau10medium = c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
                    "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
                    "#CDCC5D", "#6DCCDA")
tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
              "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
              "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
              "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5")

### =======

## Read in raw UMI x barcode matrix - **use pre-mRNA-aligned reads
path.5558hb.s3e <- file.path("/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/07_cellranger/Br5558/outs/raw_feature_bc_matrix/")
s3e.hb.5558 <- read10xCounts(path.5558hb.s3e, col.names=TRUE)
s3e.hb.5558

path.1204hb.s3e <- file.path("/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/07_cellranger/Br1204/outs/raw_feature_bc_matrix/")
s3e.hb.1204 <- read10xCounts(path.1204hb.s3e, col.names=TRUE)
s3e.hb.1204
path.5639hb.s3e <- file.path("/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/07_cellranger/Br5639/outs/raw_feature_bc_matrix/")
s3e.hb.5639 <- read10xCounts(path.5639hb.s3e, col.names=TRUE)
s3e.hb.5639

path.1735hb.s3e <- file.path("/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/07_cellranger/Br1735/outs/raw_feature_bc_matrix/")
s3e.hb.1735 <- read10xCounts(path.1735hb.s3e, col.names=TRUE)
s3e.hb.1735

path.1092hb.s3e <- file.path("/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/07_cellranger/Br1092/outs/raw_feature_bc_matrix/")
s3e.hb.1092 <- read10xCounts(path.1092hb.s3e, col.names=TRUE)
s3e.hb.1092

path.5555hb.s3e <- file.path("/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/07_cellranger/Br5555/outs/raw_feature_bc_matrix/")
s3e.hb.5555 <- read10xCounts(path.5555hb.s3e, col.names=TRUE)
s3e.hb.5555

path.1469hb.s3e <- file.path("/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/07_cellranger/Br1469/outs/raw_feature_bc_matrix/")
s3e.hb.1469 <- read10xCounts(path.1469hb.s3e, col.names=TRUE)
rm(list=ls(pattern="path."))
s3e.hb.1469


## match rownames and append by column
s3e.hb.1204 <- s3e.hb.1204[match(rownames(s3e.hb.5558), rownames(s3e.hb.1204)), ]
s3e.hb.1469 <- s3e.hb.1469[match(rownames(s3e.hb.5558), rownames(s3e.hb.1469)), ]
s3e.hb.5639 <- s3e.hb.5639[match(rownames(s3e.hb.5558), rownames(s3e.hb.5639)), ]
s3e.hb.1735 <- s3e.hb.1735[match(rownames(s3e.hb.5558), rownames(s3e.hb.1735)), ]
s3e.hb.1092 <- s3e.hb.1092[match(rownames(s3e.hb.5558), rownames(s3e.hb.1092)), ]
s3e.hb.5555 <- s3e.hb.5555[match(rownames(s3e.hb.5558), rownames(s3e.hb.5555)), ]
s3e.hb.5558 <- s3e.hb.5558[match(rownames(s3e.hb.1204), rownames(s3e.hb.5558)), ]
#######################################
###########Drop Empties###############
######################################
cat(paste0("Simulating empty drops for: s3e.hb... \n"))
set.seed(109)

e.out.hb.1204 <- emptyDrops(counts(s3e.hb.1204), niters=20000)
e.out.hb.1469 <- emptyDrops(counts(s3e.hb.1469), niters=20000)
e.out.hb.5639 <- emptyDrops(counts(s3e.hb.5639), niters=20000)
e.out.hb.1735 <- emptyDrops(counts(s3e.hb.1735), niters=20000)
e.out.hb.1092 <- emptyDrops(counts(s3e.hb.1092), niters=20000)
e.out.hb.5555 <- emptyDrops(counts(s3e.hb.5555), niters=20000)
e.out.hb.5558 <- emptyDrops(counts(s3e.hb.5558), niters=20000)
####################################
########subset for empty Drops #####
####################################

s3e.hb.1204 <- s3e.hb.1204[ ,which(e.out.hb.1204$FDR <= 0.001)]
s3e.hb.1469 <- s3e.hb.1469[ ,which(e.out.hb.1469$FDR <= 0.001)]
s3e.hb.5639 <- s3e.hb.5639[ ,which(e.out.hb.5639$FDR <= 0.001)]
s3e.hb.1735 <- s3e.hb.1735[ ,which(e.out.hb.1735$FDR <= 0.001)]
s3e.hb.1092 <- s3e.hb.1092[ ,which(e.out.hb.1092$FDR <= 0.001)]
s3e.hb.5555 <- s3e.hb.5555[ ,which(e.out.hb.5555$FDR <= 0.001)]
s3e.hb.5558 <- s3e.hb.5558[ ,which(e.out.hb.5558$FDR <= 0.001)]

sce.all.hb <- cbind(s3e.hb.1204, s3e.hb.5558, s3e.hb.1469, s3e.hb.5639, s3e.hb.1735, s3e.hb.1092, s3e.hb.5555)
rm(s3e.hb.1469,s3e.hb.1204,s3e.hb.5558,s3e.hb.5639, s3e.hb.1735, s3e.hb.1092,s3e.hb.5555)


dim(sce.all.hb )
# [1] 36601 19864

#(barcode whitelist must be ~6.8 million potential barcodes?)
rownames(sce.all.hb ) <- uniquifyFeatureNames(rowData(sce.all.hb )$ID, rowData(sce.all.hb )$Symbol)

save(sce.all.hb, e.out.hb.1204, e.out.hb.1469, e.out.hb.5639, e.out.hb.1735, e.out.hb.1092, e.out.hb.5555,
     file=here("processed-data","08_snRNA-seq_Erik","20220301_human_hb_processing.rda"))
