################################################################################
### LIBD Frankenstein S3e > 10x Chromium/library prep (> seq'g at core) test
### Initiated: MNT Nov2020
### Edited: EDN Jan2021
### Objective: To use 10x pilot pipeline and preliminarily assess data from
###            sorting with S3e & 10x Chromium/library prep, in-house
### Sample info: Br1469, habenula sorted with PI
################################################################################

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

####################################
########subset for empty Drops #####
####################################

s3e.hb.1204 <- s3e.hb.1204[ ,which(e.out.hb.1204$FDR <= 0.001)]
s3e.hb.1469 <- s3e.hb.1469[ ,which(e.out.hb.1469$FDR <= 0.001)]
s3e.hb.5639 <- s3e.hb.5639[ ,which(e.out.hb.5639$FDR <= 0.001)]
s3e.hb.1735 <- s3e.hb.1735[ ,which(e.out.hb.1735$FDR <= 0.001)]
s3e.hb.1092 <- s3e.hb.1092[ ,which(e.out.hb.1092$FDR <= 0.001)]
s3e.hb.5555 <- s3e.hb.5555[ ,which(e.out.hb.5555$FDR <= 0.001)]

s3e.hb <- cbind(s3e.hb.1204, s3e.hb.5558, s3e.hb.1469, s3e.hb.5639, s3e.hb.1735, s3e.hb.1092, s3e.hb.5555)
rm(s3e.hb.1469,s3e.hb.1204,s3e.hb.5558,s3e.hb.5639, s3e.hb.1735, s3e.hb.1092,s3e.hb.5555)


dim(s3e.hb)
# [1]   36601 1193783
#(barcode whitelist must be ~6.8 million potential barcodes?)
rownames(s3e.hb) <- uniquifyFeatureNames(rowData(s3e.hb)$ID, rowData(s3e.hb)$Symbol)

save(s3e.hb, e.out.hb.1204, e.out.hb.1469, e.out.hb.5639, e.out.hb.1735, e.out.hb.1092, e.out.hb.5555,
     file=here("processed-data","08_snRNA-seq_Erik","20220301_human_hb_processing.rda"))

### Quality control ============================================================
## - Going to ignore the adaptive NMAD-approach to outlier detection for UMI/feature count
#    because this hasn't been as straightforward in past experience (might throw away neurons)
## - Vignette for the 10x PBMC dataset (OSCA Ch.24) only does mito & droplet QC anyhow
#       - (mention that for a sample with very heterogeneous cell comp., don't want
#          to drop potential cells with low RNA content)


## Cell detection (droplet exclusion, rather)
# Can use UMI count vs barcode rank (knee/inflection plot) to decide threshold, but
#      "this unnecessarily discards libraries derived from cell types with low RNA content" (OSCA, Ch. 6)
#      -> Instead should prefer this Monte Carlo-simulation-based empty droplet test:
# Additionally:
# For any Sig==FALSE & Limited==TRUE, may need to increase n iterations (default = 10000) with 'niters='
#   - this field = whether "the computed p-value for a...barcode is bounded by the number of iterations"
cat(paste0("Simulating empty drops for: s3e.hb... \n"))
set.seed(109)
e.out.hb <- emptyDrops(counts(s3e.hb), niters=20000)

 save(s3e.hb, e.out.hb,
      file="20210524_human_hb_processing.rda")

print(table(Signif = e.out.hb$FDR <= 0.001, Limited = e.out.hb$Limited))

#Limited
#Signif  FALSE  TRUE
#FALSE 67785     0
#TRUE    718  7173
#- all are good and not lower-p-value-bound-limited

# 2,418 estimated by Cell Ranger `count`

# (and btw:)
# print(table(Signif = e.out.hb$FDR <= 0.01, Limited = e.out.hb$Limited))
#Limited
#Signif  FALSE  TRUE
#FALSE 67785     0
#TRUE    718  7173


# Subset 'nuclei called':
s3e.hb <- s3e.hb[ ,which(e.out.hb$FDR <= 0.001)]

# Save
save(s3e.hb, e.out.hb,
     file="20210525_human_hb_processing.rda")

### Mito rate QC ===
# MAD approach for mito rate thresholding
# location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(s3e.hb)$ID,
#                    column="SEQNAME", keytype="GENEID")
# ## Warning message: \n Unable to map 312 of 33538 requested IDs.
#
# ## ID those mito genes
# stats.hb <- perCellQCMetrics(s3e.hb, subsets=list(Mito=which(location=="MT")))
#
# ## Hb
# high.mito.hb <- isOutlier(stats.hb$subsets_Mito_percent, nmads=3, type="higher")
# table(high.mito.hb)
# # 2003 TRUE
# attributes(high.mito.hb) # 1.555469
# #  -> let's go with 1.25 (see comment Re: n=12 homogenate data, above)  -  373 tossed
# mitoCutoffs <- attributes(high.mito.hb)$thresholds["higher"]
# mitoCutoffs <- round(mitoCutoffs, 2)
#
#
# # Bind stats to colData
# colData(s3e.hb) <- cbind(colData(s3e.hb), stats.hb, high.mito.hb)
# colnames(colData(s3e.hb))[9] <- "high.mito"
#
# # Let's also look at low library size and low detected features
# qc.lib <- isOutlier(s3e.hb$sum, log=TRUE, type="lower")
# qc.nexprs <- isOutlier(s3e.hb$detected, log=TRUE, type="lower")
# #No outliers...so let's set a manual threshold
# qc.lib<-s3e.hb$sum < 1000
# qc.detected<-s3e.hb$detected<500
#
# discard <- qc.lib | qc.nexprs | high.mito.hb
#
#
# ##Bind library size/expressed feature logical to sce
# colData(s3e.hb) <- cbind(colData(s3e.hb), qc.lib)
#
# colData(s3e.hb) <- cbind(colData(s3e.hb), qc.nexprs)
#
# colData(s3e.hb) <- cbind(colData(s3e.hb), discard)
#
#
# ##Plot results of library size/expressed features
# pdf("/users/enelson/Habenula/plots_human/human_hb_QC_metrics_highmito_5-25-21.pdf", height=5)
# grid.arrange(
#   plotColData(s3e.hb, y="sum", colour_by="discard") +
#     scale_y_log10() + ggtitle("Total count"),
#   plotColData(s3e.hb, y="detected", colour_by="discard") +
#     scale_y_log10() + ggtitle("Detected features"),
#   plotColData(s3e.hb, y="subsets_Mito_percent",
#               colour_by="high.mito") + ggtitle(paste0("Mito % (cutoff = ", mitoCutoffs,")")),
#   ncol=3
# )
#
# dev.off()
#
# ##discard low expressed features/small libraries
# s3e.hb<-s3e.hb[,!s3e.hb$discard]
# dim(s3e.hb)
#
# save(s3e.hb, e.out.hb, high.mito.hb, file="/users/enelson/Habenula/sce_objects/human_hb_postQC_20210525.rda")

### Add in, manually, for each region (but will put here for the pipeline) === === === ===
## Aside - testing add of rowRanges data to SCE's - from [ST_project_dir]/Layer_Notebook.R
#         (see https://github.com/LieberInstitute/HumanPilot/blob/c8a3a31b991081d656ededee59da45aa0494b334/Analysis/Layer_Notebook.R#L78-L87)


## Read in the gene information from the annotation GTF file
if (verbose) message(Sys.time(), " rtracklayer::import: reading the reference GTF file")
gtf <- rtracklayer::import(reference_gtf)
gtf <- gtf[gtf$type == "gene"]
names(gtf) <- gtf$gene_id

## Match the genes
if (verbose) message(Sys.time(), " adding gene information to the SPE object")
match_genes <- match(rownames(spe), gtf$gene_id)

if (all(is.na(match_genes))) {
## Protect against scenario where one set has GENCODE IDs and the other one has ENSEMBL IDs.
warning("Gene IDs did not match. This typically happens when you are not using the same GTF file as the one that was used by SpaceRanger. For example, one file uses GENCODE IDs and the other one ENSEMBL IDs. read10xVisiumWrapper() will try to convert them to ENSEMBL IDs.", call. = FALSE)
        match_genes <- match(gsub("\\..*", "", rownames(spe)), gsub("\\..*", "", gtf$gene_id))
    }

    if (any(is.na(match_genes))) {
        warning("Dropping ", sum(is.na(match_genes)), " out of ", length(match_genes), " genes for which we don't have information on the reference GTF file. This typically happens when you are not using the same GTF file as the one that was used by SpaceRanger.", call. = FALSE)
        ## Drop the few genes for which we don't have information
        spe <- spe[!is.na(match_genes), ]
        match_genes <- match_genes[!is.na(match_genes)]
    }

    ## Keep only some columns from the gtf
    mcols(gtf) <- mcols(gtf)[, gtf_cols[gtf_cols %in% colnames(mcols(gtf))]]

    ## Add the gene info to our SPE object
    rowRanges(spe) <- gtf[match_genes]

library(rtracklayer)

## get annotation
map = read.delim("/dcl02/lieber/ajaffe/ErikNelson/Habenula/Br1204/outs/raw_feature_bc_matrix/features.tsv.gz",
                 as.is=TRUE, header=FALSE,col.names=c("EnsemblID", "Symbol", "Type"))
## get GTF, this seems like what they used
gtf = import("/dcl01/ajaffe/data/lab/singleCell/cellranger_reference/refdata-gex-GRCh38-2020-A/genes/genes.gtf")
## of length 2565061
gtf = gtf[gtf$type == "gene"]
## of length 33538
names(gtf) = gtf$gene_id
table(names(gtf) == map$EnsemblID)
#TRUE
#33538      - (so the reordering below isn't actually necessary, but safe)
gtf = gtf[map$EnsemblID]
seqlevels(gtf)[1:25] = paste0("chr", seqlevels(gtf)[1:25])
mcols(gtf) = mcols(gtf)[,c(5:9)]

table(gtf$gene_type)
##      antisense       IG_C_gene IG_C_pseudogene       IG_D_gene       IG_J_gene
#           5497              14               9              37              18
#IG_J_pseudogene       IG_V_gene IG_V_pseudogene         lincRNA  protein_coding
#              3             144             188            7484           19912
#      TR_C_gene       TR_D_gene       TR_J_gene TR_J_pseudogene       TR_V_gene
#              6               4              79               4             106
#TR_V_pseudogene
#             33

# Then to add:
table(gtf$gene_id == rowData(s3e.hb)$ID)  # all 36601 TRUE

rowData(s3e.hb) <- cbind(rowData(s3e.hb), mcols(gtf))

## Save
save(s3e.hb, e.out.hb, high.mito.hb, file="/users/enelson/Habenula/sce_objects/human_hb_postQC_20210525.rda")


### Normalization===============================================================
set.seed(13700)
clusters <- quickCluster(s3e.hb)
s3e.hb <- computeSumFactors(s3e.hb, cluster=clusters)
s3e.hb <- logNormCounts(s3e.hb)

### Dimensionality reduction/feature selection-Erik===============================================
#This is based on fitting Poisson model and then calculating Pearson residuals
#
#Let's plot highly deviant genes first; this might give us an idea about how
#many genes to retain
s3e.hb<-devianceFeatureSelection(s3e.hb, assay="counts",
                                 fam="poisson",sorted=TRUE)
plot(rowData(s3e.hb)$poisson_deviance, type="l", xlab="ranked genes",
     ylab="poisson deviance", main="Feature Selection with Deviance",ylim=c(0,100000))
hdg<-rownames(counts(s3e.hb))[1:10000]

#Run nullResiduals to calculate Pearson residuals from poisson model
s3e.hb<-nullResiduals(s3e.hb, assay = "counts",
                      fam = "poisson", type = "pearson")

##Now, PCA on residuals
set.seed(1000)
s3e.hb <- runPCA(s3e.hb,exprs_values="poisson_pearson_residuals",ncomponents=50,
                 subset_row=hdg,BSPARAM=BiocSingular::RandomParam())
reducedDimNames(s3e.hb)

set.seed(100000)
s3e.hb <- runUMAP(s3e.hb, dimred="PCA",subset_row=hdg)

#Verify length
ncol(reducedDim(s3e.hb, "PCA"))

# Save into a new data file, which will dedicate for downstream [prelim] analysis
save(s3e.hb, hdg,
     file="20210524_human_hb_downstreamAnalysis_6samples.rda")

### Clustering-Erik=============================================================
g20 <- buildSNNGraph(s3e.hb, k=20, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g20)$membership
colData(s3e.hb)$label <- factor(clust)
table(colLabels(s3e.hb))

g50 <- buildSNNGraph(s3e.hb, k=50, use.dimred = 'PCA')
clust50 <- igraph::cluster_walktrap(g50)$membership
colData(s3e.hb)$k_50_label <- factor(clust50)
table(colData(s3e.hb)$k_50_label)

g10 <- buildSNNGraph(s3e.hb, k=10, use.dimred = 'PCA')
clust10 <- igraph::cluster_walktrap(g10)$membership
colData(s3e.hb)$k_10_label <- factor(clust10)
table(colData(s3e.hb)$k_10_label)

g5 <- buildSNNGraph(s3e.hb, k=5, use.dimred = 'PCA')
clust5 <- igraph::cluster_walktrap(g5)$membership
colData(s3e.hb)$k_5_label <- factor(clust5)
table(colData(s3e.hb)$k_5_label)

markers.custom = list(
  'neuron' = c('SYT1'),# 'SNAP25', 'GRIN1','MAP2'),
  'excitatory_neuron' = c('SLC17A6'),# 'SLC17A7', 'SLC17A8'),
  'inhibitory_neuron' = c('GAD1', 'GAD2'), #'SLC32A1'),
  'mediodorsal thalamus'= c('EPHA4','PDYN', 'LYPD6B', 'LYPD6', 'S1PR1', 'GBX2', 'RAMP3', 'COX6A2', 'SLITRK6', 'DGAT2'),
  'Hb neuron specific'= c('POU2F2','POU4F1','GPR151','CALB2'),#,'GPR151','POU4F1','STMN2','CALB2','NR4A2','VAV2','LPAR1'),
  'MHB neuron specific' = c('TAC1','CHAT','CHRNB4'),#'TAC3','SLC17A7'
  'LHB neuron specific' = c('HTR2C','MMRN1'),#'RFTN1'
  'oligodendrocyte' = c('MOBP'),# 'MBP', 'PLP1'),
  'oligodendrocyte_precursor' = c('PDGFRA'),# 'VCAN', 'CSPG4', 'GPR17'),
  'microglia' = c('C3'),# 'CSF1R', 'C3'),
  'astrocyte' = c('GFAP')#,# 'AQP4'),


#  'Hb microglia' = c('CD68'),
#  'Hb astrocytes' = c('KCNJ10')
)

pdf("20210526_human_hb_markers.pdf", height=6, width=12)
for(i in 1:length(markers.custom)){
  print(
    plotExpression(s3e.hb, exprs_values = "logcounts", features=c(markers.custom[[i]]),
                   x="k_50_label", colour_by="k_50_label", point_alpha=0.5, point_size=.7,
                   add_legend=F) +
      stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar",
                   width = 0.3) +
      ggtitle(label=paste0(names(markers.custom)[i], " markers"))
  )
}
dev.off()

pdf("20210119_human_hb_markers_k10.pdf", height=6, width=12)
for(i in 1:length(markers.custom)){
  print(
    plotExpression(s3e.hb,features=c(markers.custom[[i]]),
                   x="k_10_label", colour_by="k_10_label", point_alpha=0.5, point_size=.7,add_legend=F)+
      stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar",
                   width = 0.3) +
      ggtitle(label=paste0(names(markers.custom)[i], " markers"))
  )}
dev.off()

pdf("20210119_human_hb_markers_k5.pdf", height=6, width=12)
for(i in 1:length(markers.custom)){
  print(
    plotExpression(s3e.hb,features=c(markers.custom[[i]]),
                   x="k_5_label", colour_by="k_5_label", point_alpha=0.5, point_size=.7,add_legend=F)+
      stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar",
                   width = 0.3) +
      ggtitle(label=paste0(names(markers.custom)[i], " markers"))
  )}
dev.off()

table(s3e.hb$k_5_label)
##1   2   3   4   5   6   7   8   9  10  11  12  13
#211  45 571 293 255 465  62 161   7  95  23  44  14

cellIdx <- splitit(s3e.hb$label)
sapply(cellIdx, function(x){quantile(s3e.hb[ ,x]$sum)})
#1     2       3     4       5     6        7     8      9    10
#0%     761.0  3620  1529.0  1905  1347.0  2207  1330.00   822 2590.0  7640
#25%   1988.5  5581  5764.5 11408  4321.5  8571  3621.25  1575 4298.5 11049
#50%   2919.0  6603  7441.0 14994  5864.0 10596  5156.00  2268 6459.0 12564
#75%   4388.0  8530  9229.0 18147  7578.0 13053  6818.00  3331 7891.5 15386
#100% 19202.0 19583 40036.0 39245 20181.0 51866 17324.00 22039 8920.0 58153
#11       12       13
#0%     944   852.00 17590.00
#25%   2100  1143.25 26993.50
#50%   3265  1587.50 30555.00
#75%   4180  1942.50 34731.25
#100% 23960 22807.00 53094.00

##make level 1
colData(s3e.hb)$level_1<-factor(
  ifelse(colData(s3e.hb)$label %in% c(2,7,9),'LHb_neuron',
      ifelse(colData(s3e.hb)$label %in% c(12,13),'MHb_neuron',
            ifelse(colData(s3e.hb)$label %in% c(3,10),'Oligodendrocyte',
                   ifelse(colData(s3e.hb)$label %in% c(5,15),'Ambiguous',
                        ifelse(colData(s3e.hb)$label %in% c(4,14),'Thalamus_neuron',
                               ifelse(colData(s3e.hb)$label %in% c(8,11),'OPC',
                                      ifelse(colData(s3e.hb)$label==1,'Microglia',
                                             'Astrocyte')
)))))))

#make sample_name
colData(s3e.hb)$sample_name<-factor(
  ifelse(colData(s3e.hb)$Sample == "/dcl02/lieber/ajaffe/ErikNelson/Habenula/Br1204/outs/raw_feature_bc_matrix/",'Br1204',
         ifelse(colData(s3e.hb)$Sample == "/dcl02/lieber/ajaffe/ErikNelson/Habenula/Br1469/outs/raw_feature_bc_matrix/",'Br1469',
                'Br5558')))
save(s3e.hb,hdg,marks,marks.lvl1,file="human_hb_downstreamAnalysis_3samples.rda")

##Plot heatmap of interesting genes
s3e.hb<-s3e.hb[,unique(colnames(s3e.hb))]
nicotine.genes<-c("CHRNA3", "CHRNA5", "CHRNB4", "CHRNA2", "CHRNB2", "TCF7L2", "GPR151",
         "DRD3", "GLP1R")
opioid.genes<-c("CHRNA3", "CHRNA5", "CHRNB4","OPRM1", "TGFBR1", "TLR4", "GPR139", "RPS6KA3","CAMK2A","CAMK2B")
alcohol.genes<-c("CHRNA3", "CHRNA5", "CHRNB4","GLRA1", "GLRA2", "GLRA3", "GLRA4", "GLRB", "GRM5", "TRPV1", "HTR2C", "KCNQ2", "KCNQ3","CAMK2A","CAMK2B")
cocaine.genes<-c("CHRNA3", "CHRNA5", "CHRNB4","NR4A2")

pdf("nicotine_genes_heatmap.pdf")
plotHeatmap(s3e.subset,features=nicotine.genes,order_columns_by="label",
            colour_columns_by="label",main="Nicotine-associated gene expression")
dev.off()

pdf("opioid_genes_heatmap.pdf")
plotHeatmap(s3e.subset,features=opioid.genes,order_columns_by="label",
            colour_columns_by="label", main="Opioid-associated gene expression")
dev.off()

pdf("alcohol_genes_heatmap.pdf")
plotHeatmap(s3e.subset,features=alcohol.genes,order_columns_by="label",
            colour_columns_by="label", main="Alcohol-associated gene expression")
dev.off()

pdf("cocaine_genes_heatmap.pdf")
plotHeatmap(s3e.subset,features=cocaine.genes,order_columns_by="label",
            colour_columns_by="label", main="Cocaine-associated gene expression")
dev.off()

colData(s3e.hb)$cell.type.broad<-factor(
  ifelse(s3e.hb$label==1, "LHb.neuron.1",
         ifelse(s3e.hb$label==6, "LHb.neuron.2",
                ifelse(s3e.hb$label==16, "LHb.neuron.3",
                       ifelse(s3e.hb$label==5, "LHb.neuron.4",
                              ifelse(s3e.hb$label==10, "MHb.neuron.1",
                                     ifelse(s3e.hb$label==14, "MHb.neuron.2",
                                            ifelse(s3e.hb$label==2, "Thal.neuron.1",
                                                   ifelse(s3e.hb$label==7, "Thal.neuron.2",
                                                          ifelse(s3e.hb$label==15, "Thal.neuron.3",
                                                                 ifelse(s3e.hb$label==3, "Astrocyte",
                                                                        ifelse(s3e.hb$label==4, "Microglia",
                                                                               ifelse(s3e.hb$label==9, "OPC.1",
                                                                                      ifelse(s3e.hb$label==13, "OPC.2",
                                                                                             ifelse(s3e.hb$label==8, "Oligo.1",
                                                                                                    ifelse(s3e.hb$label==11, "Oligo.2",
                                                                                                           "Oligo.3"
                                                                                                    ))))))))))))))))
## MNT wrap up prelim analysis 20Nov2020 =======================================================

## Plot PCA
pdf("pdfs/zS3e-v2-assessment_PCA_Hb-and-AMY-NeuN_prelimClusters_Nov2020.pdf", height=7, width=7)
# AMY-NeuN
plotReducedDim(s3e.amy.neun, dimred="PCA", ncomponents=4, colour_by="prelimCluster", point_alpha=0.5) +
  ggtitle(label="PI,NeuN-sorted AMY: Top 4 PCs")
# Hb
plotReducedDim(s3e.hb, dimred="PCA", ncomponents=4, colour_by="prelimCluster", text_by="prelimCluster",
               point_alpha=0.5) + ggtitle(label="PI-sorted Hb: Top 4 PCs")
dev.off()


### Generate some cluster-vs-all t's, using traditional t-test in `findMarkers()`=== ===

## Hb ===

# Remove 0 genes across all nuclei
s3e.hb <- s3e.hb[!rowSums(assay(s3e.hb, "counts"))==0, ]  # keeps 25904 genes

#mod <- with(colData(sce.hb.st), model.matrix(~ donor))      - not relevant; only one sample

markers.hb.s3e <- list()
for(i in levels(s3e.hb$prelimCluster)){
  # Make temporary contrast
  s3e.hb$contrast <- ifelse(s3e.hb$prelimCluster==i, 1, 0)
  # Test cluster vs. all
  markers.hb.s3e[[i]] <- findMarkers(s3e.hb, groups=s3e.hb$contrast,
                                     assay.type="logcounts", test="t", #design=mod,
                                     direction="up", pval.type="all", full.stats=T)
}

## Temp set of stats to get the standardized logFC
temp.1vAll <- list()
for(i in levels(s3e.hb$prelimCluster)){
  # Make temporary contrast
  s3e.hb$contrast <- ifelse(s3e.hb$prelimCluster==i, 1, 0)
  # Test cluster vs. all
  temp.1vAll[[i]] <- findMarkers(s3e.hb, groups=s3e.hb$contrast,
                                 assay.type="logcounts", test="t", #design=mod,
                                 std.lfc=TRUE,
                                 direction="up", pval.type="all", full.stats=T)
}


# Replace that empty slot with the entry with the actul stats
markers.hb.s3e <- lapply(markers.hb.s3e, function(x){ x[[2]] })
# Same for that with std.lfc
temp.1vAll <- lapply(temp.1vAll, function(x){ x[[2]] })

# Now just pull from the 'stats.0' DataFrame column
markers.hb.s3e <- lapply(markers.hb.s3e, function(x){ x$stats.0 })
temp.1vAll <- lapply(temp.1vAll, function(x){ x$stats.0 })

# Re-name std.lfc column and add to the first result
for(i in names(temp.1vAll)){
  colnames(temp.1vAll[[i]])[1] <- "std.logFC"
  markers.hb.s3e[[i]] <- cbind(markers.hb.s3e[[i]], temp.1vAll[[i]]$std.logFC)
  # Oh the colname is kept weird
  colnames(markers.hb.s3e[[i]])[4] <- "std.logFC"
  # Then re-organize
  markers.hb.s3e[[i]] <- markers.hb.s3e[[i]][ ,c("logFC","std.logFC","log.p.value","log.FDR")]
}



### Same for AMY-NeuN sample === === ============================================
# Remove 0 genes across all nuclei
s3e.amy.neun <- s3e.amy.neun[!rowSums(assay(s3e.amy.neun, "counts"))==0, ]  # keeps 25610 genes

markers.amy.neun.s3e <- list()
for(i in levels(s3e.amy.neun$prelimCluster)){
  # Make temporary contrast
  s3e.amy.neun$contrast <- ifelse(s3e.amy.neun$prelimCluster==i, 1, 0)
  # Test cluster vs. all
  markers.amy.neun.s3e[[i]] <- findMarkers(s3e.amy.neun, groups=s3e.amy.neun$contrast,
                                           assay.type="logcounts", test="t", #design=mod,
                                           direction="up", pval.type="all", full.stats=T)
}

## Temp set of stats to get the standardized logFC
temp.1vAll <- list()
for(i in levels(s3e.amy.neun$prelimCluster)){
  # Make temporary contrast
  s3e.amy.neun$contrast <- ifelse(s3e.amy.neun$prelimCluster==i, 1, 0)
  # Test cluster vs. all
  temp.1vAll[[i]] <- findMarkers(s3e.amy.neun, groups=s3e.amy.neun$contrast,
                                 assay.type="logcounts", test="t", #design=mod,
                                 std.lfc=TRUE,
                                 direction="up", pval.type="all", full.stats=T)
}


# Replace that empty slot with the entry with the actul stats
markers.amy.neun.s3e <- lapply(markers.amy.neun.s3e, function(x){ x[[2]] })
# Same for that with std.lfc
temp.1vAll <- lapply(temp.1vAll, function(x){ x[[2]] })

# Now just pull from the 'stats.0' DataFrame column
markers.amy.neun.s3e <- lapply(markers.amy.neun.s3e, function(x){ x$stats.0 })
temp.1vAll <- lapply(temp.1vAll, function(x){ x$stats.0 })

# Re-name std.lfc column and add to the first result
for(i in names(temp.1vAll)){
  colnames(temp.1vAll[[i]])[1] <- "std.logFC"
  markers.amy.neun.s3e[[i]] <- cbind(markers.amy.neun.s3e[[i]], temp.1vAll[[i]]$std.logFC)
  # Oh the colname is kept weird
  colnames(markers.amy.neun.s3e[[i]])[4] <- "std.logFC"
  # Then re-organize
  markers.amy.neun.s3e[[i]] <- markers.amy.neun.s3e[[i]][ ,c("logFC","std.logFC","log.p.value","log.FDR")]
}


## Let's save this along with the previous pairwise results
save(markers.hb.s3e, markers.amy.neun.s3e,
     file="rdas/zS3e-v2-pilot_markers-stats_Hb-and-AMY-NeuN_findMarkers-SN-LEVEL_Nov2020.rda")





### Compare to 10x pilot datasets ========================================
library(pheatmap)
library(RColorBrewer)

## Set up S3e Hb t's
fixTo <- rownames(markers.hb.s3e[[1]])

for(s in names(markers.hb.s3e)){
  markers.hb.s3e[[s]]$t.stat <- markers.hb.s3e[[s]]$std.logFC * sqrt(ncol(s3e.hb))
  markers.hb.s3e[[s]] <- markers.hb.s3e[[s]][fixTo, ]
}

# Pull out the t's
ts.hb.s3e <- sapply(markers.hb.s3e, function(x){x$t.stat})
rownames(ts.hb.s3e) <- fixTo


## Pull in t's from across the brain ===
load("rdas/markers-stats_SN-LEVEL_1vAll_all-regions-combined_May2020.rda", verbose=T)
# FMstats.list, sampleNumNuclei, readme.mnt, ref.sampleInfo
readme.mnt

## Create matrix of t's with region:subcluster identifiers
ts.list <- lapply(FMstats.list, function(x){
  sapply(x, function(y){y$t.stat})
}
)
# Add back in row names and region suffix
for(i in names(ts.list)){
  rownames(ts.list[[i]]) <- rownames(FMstats.list[[1]][[1]])
  colnames(ts.list[[i]]) <- paste0(colnames(ts.list[[i]]), "_", i)
}
# Cbind
ts.fullMat <- do.call(cbind, ts.list)


# Shorten names
colnames(ts.fullMat) <- gsub("Excit", "Ex", colnames(ts.fullMat))
colnames(ts.fullMat) <- gsub("Inhib", "In", colnames(ts.fullMat))

# Intersecting expressed genes
sharedGenes.all <- intersect(rownames(ts.hb.s3e), rownames(ts.fullMat))
# of length 25,137

# Subset/order
ts.hb.s3e <- ts.hb.s3e[sharedGenes.all, ]
ts.fullMat <- ts.fullMat[sharedGenes.all, ]

colnames(ts.hb.s3e) <- paste0("clust_", colnames(ts.hb.s3e))

cor_t_hb <- cor(ts.fullMat, ts.hb.s3e)
range(cor_t_hb)


## Heatmap
theSeq.all = seq(-.90, .90, by = 0.01)
my.col.all <- colorRampPalette(brewer.pal(7, "BrBG"))(length(theSeq.all)-1)

pdf("pdfs/zS3e-v2-assessment_overlap-S3e-Hb_to_LIBD-10x-pilot-68-subclusts_Nov2020.pdf", width=6, height=10)
# By region (automatically ordered)
pheatmap(cor_t_hb,
         color=my.col.all,
         cluster_cols=F, cluster_rows=F,
         breaks=theSeq.all,
         fontsize=8, fontsize_row=10, fontsize_col=11,
         display_numbers=T, number_format="%.2f", fontsize_number=4.8,
         legend_breaks=c(seq(-0.90,0.90,by=0.45)),
         main="Correlation of cluster-specific t's: S3e (PI)-sorted Hb vs. \n 68 subclusters from (Tran et al. 2020), Nov2020")
# Ordering by name (and ~ broad cell type)
pheatmap(cor_t_hb[sort(rownames(cor_t_hb)), ],
         color=my.col.all,
         cluster_cols=F, cluster_rows=F,
         breaks=theSeq.all,
         fontsize=8, fontsize_row=10, fontsize_col=11,
         display_numbers=T, number_format="%.2f", fontsize_number=4.8,
         legend_breaks=c(seq(-0.90,0.90,by=0.45)),
         main="Correlation of cluster-specific t's: S3e (PI)-sorted Hb vs. \n 68 subclusters from (Tran et al. 2020), Nov2020")
dev.off()





## Similarly, for NAc-NeuN sample === === ===

## Set up S3e NAc sample t's
fixTo <- rownames(markers.amy.neun.s3e[[1]])

for(s in names(markers.amy.neun.s3e)){
  markers.amy.neun.s3e[[s]]$t.stat <- markers.amy.neun.s3e[[s]]$std.logFC * sqrt(ncol(s3e.amy.neun))
  markers.amy.neun.s3e[[s]] <- markers.amy.neun.s3e[[s]][fixTo, ]
}

# Pull out the t's
ts.amy.neun.s3e <- sapply(markers.amy.neun.s3e, function(x){x$t.stat})
rownames(ts.amy.neun.s3e) <- fixTo


# Intersecting expressed genes (with 'ts.fullMat', as above)
ts.fullMat <- do.call(cbind, ts.list)
sharedGenes.all <- intersect(rownames(ts.amy.neun.s3e), rownames(ts.fullMat))
# of length 24,854

# Subset/order
ts.amy.neun.s3e <- ts.amy.neun.s3e[sharedGenes.all, ]
ts.fullMat <- ts.fullMat[sharedGenes.all, ]

colnames(ts.amy.neun.s3e) <- paste0("clust_", colnames(ts.amy.neun.s3e))

cor_t_amy <- cor(ts.fullMat, ts.amy.neun.s3e)
range(cor_t_amy)


## Heatmap
theSeq.all = seq(-.90, .90, by = 0.01)
my.col.all <- colorRampPalette(brewer.pal(7, "BrBG"))(length(theSeq.all)-1)

pdf("pdfs/zS3e-v2-assessment_overlap-S3e-AMY-NeuN_to_LIBD-10x-pilot-68-subclusts_Nov2020.pdf", width=7, height=10)
# By region (automatically ordered)
pheatmap(cor_t_amy,
         color=my.col.all,
         cluster_cols=F, cluster_rows=F,
         breaks=theSeq.all,
         fontsize=8, fontsize_row=10, fontsize_col=11,
         display_numbers=T, number_format="%.2f", fontsize_number=4.8,
         legend_breaks=c(seq(-0.90,0.90,by=0.45)),
         main="Correlation of cluster-specific t's: S3e (PI/NeuN)-sorted Hb vs. \n 68 subclusters from (Tran et al. 2020), Nov2020")
# Ordering by name (and ~ broad cell type)
pheatmap(cor_t_amy[sort(rownames(cor_t_amy)), ],
         color=my.col.all,
         cluster_cols=F, cluster_rows=F,
         breaks=theSeq.all,
         fontsize=8, fontsize_row=10, fontsize_col=11,
         display_numbers=T, number_format="%.2f", fontsize_number=4.8,
         legend_breaks=c(seq(-0.90,0.90,by=0.45)),
         main="Correlation of cluster-specific t's: S3e (PI/NeuN)-sorted Hb vs. \n 68 subclusters from (Tran et al. 2020), Nov2020")
dev.off()

##do fastMNN() on both new SCE objects
s3e.hb$Sample<-factor(s3e.hb$Sample)
set.seed(1000101001)
mnn.out <- fastMNN(s3e.hb, batch=s3e.hb$Sample, d=50, k=20, subset.row=hdg,
                   assay.type="logcounts",
                   BSPARAM=BiocSingular::RandomParam(deferred=TRUE))
mnn.out

set.seed(100000)
mnn.out <- runUMAP(mnn.out, dimred="corrected")

g20 <- buildSNNGraph(mnn.out, k=20, use.dimred = 'corrected')
clust <- igraph::cluster_walktrap(g20)$membership
colData(mnn.out)$label <- factor(clust)
table(colLabels(mnn.out))

g50 <- buildSNNGraph(mnn.out, k=50, use.dimred = 'corrected')
clust50 <- igraph::cluster_walktrap(g50)$membership
colData(mnn.out)$k_50_label <- factor(clust50)
table(colData(mnn.out)$k_50_label)

plotUMAP(mnn.out,text_by='label',colour_by='batch',by_exprs_features='reconstructed')
plotUMAP(mnn.out,text_by='label',colour_by='label',by_exprs_features='reconstructed')
plotUMAP(mnn.out,text_by='label',colour_by='SYT1',by_exprs_features='reconstructed')
plotUMAP(mnn.out,text_by='label',colour_by='LYPD6B',by_exprs_features='reconstructed')
plotUMAP(mnn.out,text_by='label',colour_by='HTR2C',by_exprs_features='reconstructed')
plotUMAP(mnn.out,text_by='label',colour_by='CHAT',by_exprs_features='reconstructed')
plotUMAP(mnn.out,text_by='label',colour_by='TAC1',by_exprs_features='reconstructed')
plotUMAP(mnn.out,text_by='k_50_label',colour_by='k_50_label',by_exprs_features='reconstructed')


features<-c("GFAP","CSF1R","MOBP","PDGFRA","SNAP25","GPR151","HTR2C","CHRNB3","LYPD6B","ADARB2","NECAB1")
plotHeatmap(s3e.hb,features=features,order_columns_by="cell.type",
            colour_columns_by="cell.type",cluster_rows=F)

plotHeatmap(s3e.hb,features=features,order_columns_by="cell.type",
            colour_columns_by="cell.type",cluster_rows=F)

plotGroupedHeatmap(s3e.hb,features=features,cluster_rows=F,cluster_cols=F,group='cell.type')

g50 <- buildSNNGraph(s3e.hb, k=50, use.dimred = 'corrected')
clust50 <- igraph::cluster_walktrap(g50)$membership
colData(s3e.hb)$k_50_label <- factor(clust50)
table(colData(s3e.hb)$k_50_label)

'PDYN', 'LYPD6B', 'LYPD6', 'S1PR1', 'GBX2', 'RAMP3', 'COX6A2', 'SLITRK6', 'DGAT2'

colData(s3e.hb)$cell.type<-factor(
  ifelse(s3e.hb$k_50_label %in% c(16,18,25), 'Astrocyte',
         ifelse(s3e.hb$k_50_label %in% c(4,13), 'Microglia',
                ifelse(s3e.hb$k_50_label %in% c(6,17,23), 'OPC',
                       ifelse(s3e.hb$k_50_label %in% c(1,2,24), 'Oligodendrocyte',
                              ifelse(s3e.hb$k_50_label %in% c(12,20), 'Neuron.MHb',
                                     ifelse(s3e.hb$k_50_label %in% c(7,8,9,10,22), 'Neuron.Thalamus.MD',
                                            ifelse(s3e.hb$k_50_label %in% c(15,21,26),'Neuron.Thalamus.PV/PF',
                                                   'Neuron.LHb'))))))))
