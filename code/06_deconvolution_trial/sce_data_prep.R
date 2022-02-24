library(SingleCellExperiment)
library(jaffelab)
library(here)
library(GenomicFeatures)
library(sessioninfo)

#### Load and filter data ####
## Load rse_gene data
set_here(path='..')
load(here("count_data","rse_gene_Roche_Habenula_qcAndAnnotated_n69.Rdata"), verbose = TRUE)
rownames(rse_gene) <- rowData(rse_gene)$ensemblID

## All data
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/all-FACS-homogenates_n12_PAN-BRAIN-Analyses_MNTFeb2020.rda", verbose = TRUE)
sce.all.n12$uniqueID <- paste0(sce.all.n12$donor, "_", sce.all.n12$Barcode)
colnames(sce.all.n12) <- sce.all.n12$uniqueID

sce.all.n12 <- sce.all.n12[,sce.all.n12$cellType != "Ambig.lowNtrxts",]
sce.all.n12$cellType <- droplevels(sce.all.n12$cellType)
## Add cellType.broad
sce.all.n12$cellType.Broad <- ss(as.character(sce.all.n12$cellType), "\\.", 1)
sce.all.n12$cellType.Broad <- factor(sce.all.n12$cellType.Broad)
## Match rownames
rownames(sce.all.n12) <- rowData(sce.all.n12)$ID
table(rownames(sce.all.n12) %in% rownames(rse_gene))
#FALSE  TRUE
#  341 33197

## filter for common between sc and bulk data
common_genes <- rownames(rse_gene)[rownames(rse_gene) %in% rownames(sce.all.n12)]
sce.all.n12 <- sce.all.n12[common_genes, ]

# Remove 0 genes across all nuclei
sce.all.n12 <- sce.all.n12[!rowSums(assay(sce.all.n12, "counts"))==0, ]
dim(sce.all.n12)
# [1] 29818 34070

## Amyg Data

#### Add bp length for RPKM later ####
rd.all <- rowData(sce.all.n12)
## Import gnomic features
txdb <- makeTxDbFromGFF("/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/GRCh38-3.0.0_premrna/genes/genes.gtf")
g <- genes(txdb)
e <- exonsBy(txdb, by = "gene")

if (!all(names(e) %in% names(g))) {
  warning("Dropping exons with gene ids not present in the gene list")
  e <- e[names(e) %in% names(g)]
}
e2 <- disjoin(e)
g$bp_length <- sum(width(e2))
summary(g$bp_length)

## Add data to sacc
g_all <- g[rownames(sce.all.n12),]
summary(g_all$bp_length)
table(g_all$bp_length > 205012)

rowRanges(sce.all.n12) <- g_all
all(rownames(rd.all) == rownames(sce.all.n12))
rowData(sce.all.n12)$Symbol <- rd.all$Symbol


## Save filtered sce object
save(sce.all.n12, file = here("Deconvolution","data","sce.all.n12_filtered.Rdata"))



# sgejobs::job_single('sce_data_prep', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript 1_sce_data_prep.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# [1] "Reproducibility information:"
# > Sys.time()
# [1] "2021-01-21 13:16:34 EST"
# > proc.time()
#     user   system  elapsed
#  307.517   18.237 2664.127
# > options(width = 120)
# > session_info()
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.0.2 Patched (2020-06-24 r78746)
#  os       CentOS Linux 7 (Core)
# 0;38;5;40m CRAN (R 4.0.2)
#  Biobase              * 2.48.0   2020-04-27 [2] Bioconductor
#  BiocFileCache          1.12.1   2020-08-04 [2] Bioconductor
#  BiocGenerics         * 0.34.0   2020-04-27 [2] Bioconductor
#  BiocParallel           1.22.0   2020-04-27 [2] Bioconductor
#  biomaRt                2.44.1   2020-06-17 [1] Bioconductor
#  Biostrings             2.56.0   2020-04-27 [2] Bioconductor
#  bit                    4.0.4    2020-08-04 [2] CRAN (R 4.0.2)
#  bit64                  4.0.5    2020-08-30 [2] CRAN (R 4.0.2)
#  bitops                 1.0-6    2013-08-17 [2] CRAN (R 4.0.0)
#  blob                   1.2.1    2020-01-20 [2] CRAN (R 4.0.0)
#  cli                    2.1.0    2020-10-12 [2] CRAN (R 4.0.2)
#  colorout             * 1.2-2    2020-06-01 [1] Github (jalvesaq/colorout@726d681)
#  crayon                 1.3.4    2017-09-16 [2] CRAN (R 4.0.0)
#  curl                   4.3      2019-12-02 [2] CRAN (R 4.0.0)
#  DBI                    1.1.0    2019-12-15 [2] CRAN (R 4.0.0)
#  dbplyr                 1.4.4    2020-05-27 [2] CRAN (R 4.0.2)
#  DelayedArray         * 0.14.1   2020-07-14 [2] Bioconductor
#  digest                 0.6.26   2020-10-17 [2] CRAN (R 4.0.2)
#  dplyr                  1.0.2    2020-08-18 [2] CRAN (R 4.0.2)
#  ellipsis               0.3.1    2020-05-15 [2] CRAN (R 4.0.0)
#  fansi                  0.4.1    2020-01-08 [2] CRAN (R 4.0.0)
#  generics               0.0.2    2018-11-29 [2] CRAN (R 4.0.0)
#  GenomeInfoDb         * 1.24.2   2020-06-15 [2] Bioconductor
#  GenomeInfoDbData       1.2.3    2020-05-18 [2] Bioconductor
#  GenomicAlignments      1.24.0   2020-04-27 [2] Bioconductor
#  GenomicFeatures      * 1.40.1   2020-07-08 [2] Bioconductor
#  GenomicRanges        * 1.40.0   2020-04-27 [2] Bioconductor
#  glue                   1.4.2    2020-08-27 [2] CRAN (R 4.0.2)
#  googledrive            1.0.1    2020-05-05 [1] CRAN (R 4.0.0)
#  here                 * 0.1-11   2020-08-07 [1] Github (krlmlr/here@d0feb09)
#  hms                    0.5.3    2020-01-08 [2] CRAN (R 4.0.0)
#  httr                   1.4.2    2020-07-20 [2] CRAN (R 4.0.2)
#  IRanges              * 2.22.2   2020-05-21 [2] Bioconductor
#  jaffelab             * 0.99.30  2020-06-03 [1] Github (LieberInstitute/jaffelab@42637ff)
#  lattice                0.20-41  2020-04-02 [3] CRAN (R 4.0.2)
#  lifecycle              0.2.0    2020-03-06 [2] CRAN (R 4.0.0)
#  limma                  3.44.3   2020-06-12 [2] Bioconductor
#  magrittr               1.5      2014-11-22 [2] CRAN (R 4.0.0)
#  Matrix                 1.2-18   2019-11-27 [3] CRAN (R 4.0.2)
#  matrixStats          * 0.57.0   2020-09-25 [2] CRAN (R 4.0.2)
#  memoise                1.1.0    2017-04-21 [2] CRAN (R 4.0.0)
#  openssl                1.4.3    2020-09-18 [2] CRAN (R 4.0.2)
#  pillar                 1.4.6    2020-07-10 [2] CRAN (R 4.0.2)
#  pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.0.0)
#  prettyunits            1.1.1    2020-01-24 [2] CRAN (R 4.0.0)
#  progress               1.2.2    2019-05-16 [2] CRAN (R 4.0.0)
#  purrr                  0.3.4    2020-04-17 [2] CRAN (R 4.0.0)
#  R6                     2.4.1    2019-11-12 [2] CRAN (R 4.0.0)
#  rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 4.0.0)
#  rappdirs               0.3.1    2016-03-28 [2] CRAN (R 4.0.0)
#  RColorBrewer           1.1-2    2014-12-07 [2] CRAN (R 4.0.0)
#  Rcpp                   1.0.5    2020-07-06 [2] CRAN (R 4.0.2)
#  RCurl                  1.98-1.2 2020-04-18 [2] CRAN (R 4.0.0)
#  rlang                  0.4.8    2020-10-08 [2] CRAN (R 4.0.2)
#  rprojroot              1.3-2    2018-01-03 [2] CRAN (R 4.0.0)
#  Rsamtools              2.4.0    2020-04-27 [2] Bioconductor
#  RSQLite                2.2.1    2020-09-30 [2] CRAN (R 4.0.2)
#  rtracklayer            1.48.0   2020-04-27 [2] Bioconductor
#  S4Vectors            * 0.26.1   2020-05-16 [2] Bioconductor
#  segmented              1.2-0    2020-06-23 [2] CRAN (R 4.0.2)
#  sessioninfo          * 1.1.1    2018-11-05 [2] CRAN (R 4.0.0)
#  SingleCellExperiment * 1.10.1   2020-04-28 [2] Bioconductor
#  stringi                1.5.3    2020-09-09 [2] CRAN (R 4.0.2)
#  stringr                1.4.0    2019-02-10 [2] CRAN (R 4.0.0)
#  SummarizedExperiment * 1.18.2   2020-07-09 [2] Bioconductor
#  tibble                 3.0.4    2020-10-12 [2] CRAN (R 4.0.2)
#  tidyselect             1.1.0    2020-05-11 [2] CRAN (R 4.0.0)
#  vctrs                  0.3.4    2020-08-29 [2] CRAN (R 4.0.2)
#  withr                  2.3.0    2020-09-22 [2] CRAN (R 4.0.2)
#  XML                    3.99-0.5 2020-07-23 [2] CRAN (R 4.0.2)
#  XVector                0.28.0   2020-04-27 [2] Bioconductor
#  zlibbioc               1.34.0   2020-04-27 [2] Bioconductor
#
# [1] /users/jstolz/R/4.0
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0/R/4.0/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0/R/4.0/lib64/R/library
#
