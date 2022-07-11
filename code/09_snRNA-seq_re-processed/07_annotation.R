library(SingleCellExperiment)
library(scRNAseq)
library(batchelor)
library(scater)
library(scran)
library(uwot)
library(DropletUtils)
library(jaffelab)
library(Rtsne)
library(here)
library(utils)
library(devtools)


load(here("processed-data","09_snRNA-seq_re-processed","05_collapsedClustering.Rda"))

n_clusters <- length(levels(sce.all.hb$collapsedCluster))

## Add annotations, looking at marker gene expression
annotationTab.hab <- data.frame(collapsedCluster=c(1:n_clusters))
annotationTab.hab$cellType <- NA
annotationTab.hab$cellType[c(4)] <- c("Inhib")
annotationTab.hab$cellType[c(1,2,5,6)] <- paste0("Excit_", c("A","B","C","D"))
annotationTab.hab$cellType[c(9,10)] <- paste0("Astro_",c("A","B"))
annotationTab.hab$cellType[c(11)] <- c("Micro")
annotationTab.hab$cellType[c(7)] <- c("OPC")
annotationTab.hab$cellType[c(3)] <- c("Oligo")
annotationTab.hab$cellType[c(8)] <- c("Ambig")

sce.all.hb$cellType <- annotationTab.hab$cellType[match(sce.all.hb$collapsedCluster,
                                                         annotationTab.hab$collapsedCluster)]

sce.all.hb$cellType <- factor(sce.all.hb$cellType)

table(sce.all.hb$cellType)
# Ambig Astro_A Astro_B Excit_A Excit_B Excit_C   Inhib   Micro   Oligo     OPC
#     225     302     183    1989    2595     692    6714     211    2168     110

n_Prelimclusters <- length(levels(sce.all.hb$prelimCluster))

## Add annotations, looking at marker gene expression
annotationTab.hab <- data.frame(prelimCluster=c(1:n_Prelimclusters))
annotationTab.hab$Region <- NA
annotationTab.hab$Region[c(9,12)] <- c("Habenula_broad")
annotationTab.hab$Region[c(11,22)] <- paste0("Medial-Hab",c("A","B"))
annotationTab.hab$Region[c(6,4,8,13,16,19)] <- paste0("Lateral-Hab_", c("A","B","C","D","E","F"))
annotationTab.hab$Region[c(15,18,20,21,23)] <- paste0("Thalmus_Broad",c("A","B","C","D","E"))
annotationTab.hab$Region[c(14)] <- c("Thalmus_PF")
annotationTab.hab$Region[c(3)] <- c("Thalmus_PVT")
annotationTab.hab$Region[c(1,2,5,7,5,10,17,24,25,26)] <- paste0("Ambig",c("A","B","C","D","E","F","G","H","I","J"))


sce.all.hb$Region <- annotationTab.hab$Region[match(sce.all.hb$prelimCluster,
                                                         annotationTab.hab$prelimCluster)]

sce.all.hb$Region <- factor(sce.all.hb$Region)

table(sce.all.hb$Region, sce.all.hb$cellType)
  #                  Ambig Astro_A Astro_B Excit_A Excit_B Excit_C Excit_D Inhib Micro Oligo  OPC
  # AmbigA           225       0       0       0       0       0       0     0     0     0    0
  # AmbigB             0       0       0       0       0       0       0     0     0  1115    0
  # AmbigD             0       0       0       0       0       0       0     0     0   398    0
  # AmbigE             0     302       0       0       0       0       0     0     0     0    0
  # AmbigF             0       0       0       0       0       0       0     0     0   306    0
  # AmbigG             0       0       0       0       0       0       0     0     0     0  798
  # AmbigH             0       0       0       0       0       0       0     0     0     0   79
  # AmbigI             0       0     183       0       0       0       0     0     0     0    0
  # AmbigJ             0       0       0       0       0       0       0     0     0   349    0
  # Habenula_broad     0       0       0       0    1553       0       0     0   211     0    0
  # Lateral-Hab_A      0       0       0       0     279       0       0     0     0     0    0
  # Lateral-Hab_B      0       0       0     375       0       0       0     0     0     0    0
  # Lateral-Hab_C      0       0       0     462       0       0       0     0     0     0    0
  # Lateral-Hab_D      0       0       0     899       0       0       0     0     0     0    0
  # Lateral-Hab_E      0       0       0     141       0       0       0     0     0     0    0
  # Lateral-Hab_F      0       0       0     112       0       0       0     0     0     0    0
  # Medial-HabA        0       0       0       0       0     544       0     0     0     0    0
  # Medial-HabB        0       0       0       0       0     148       0     0     0     0    0
  # Thalmus_BroadA     0       0       0       0     186       0       0     0     0     0    0
  # Thalmus_BroadB     0       0       0       0       0       0       0  3986     0     0    0
  # Thalmus_BroadC     0       0       0       0       0       0      71     0     0     0    0
  # Thalmus_BroadD     0       0       0       0       0       0      39     0     0     0    0
  # Thalmus_BroadE     0       0       0       0       0       0       0  1982     0     0    0
  # Thalmus_PF         0       0       0       0       0       0       0   746     0     0    0
  # Thalmus_PVT        0       0       0       0     577       0       0     0     0     0    0



save(sce.all.hb, file = here("processed-data","09_snRNA-seq_re-processed","07_annotation.Rda"))


### Make plots ###


pdf(here("plots","09_snRNA-seq_re-processed","regionSpecific_HAB-n7_reducedDims-with-annotated.pdf"))
plotTSNE(sce.all.hb, colour_by="Region", point_alpha=0.5)
plotUMAP(sce.all.hb, colour_by="Region", point_alpha=0.5)
plotTSNE(sce.all.hb, colour_by="cellType", point_alpha=0.5)
plotUMAP(sce.all.hb, colour_by="cellType", point_alpha=0.5)
dev.off()

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
 # "2022-07-08 14:18:59 EDT"
 #    user   system  elapsed
 # 121.953    9.835 2891.346
 # 2021-10-26 [2] Bioconductor
 # AnnotationHub            3.2.2     2022-03-01 [2] Bioconductor
 # assertthat               0.2.1     2019-03-21 [2] CRAN (R 4.1.0)
 # batchelor              * 1.10.0    2021-10-26 [1] Bioconductor
 # beachmat                 2.10.0    2021-10-26 [2] Bioconductor
 # beeswarm                 0.4.0     2021-06-01 [2] CRAN (R 4.1.2)
 # Biobase                * 2.54.0    2021-10-26 [2] Bioconductor
 # BiocFileCache            2.2.1     2022-01-23 [2] Bioconductor
 # BiocGenerics           * 0.40.0    2021-10-26 [2] Bioconductor
 # BiocIO                   1.4.0     2021-10-26 [2] Bioconductor
 # BiocManager              1.30.18   2022-05-18 [2] CRAN (R 4.1.2)
 # BiocNeighbors            1.12.0    2021-10-26 [2] Bioconductor
 # BiocParallel             1.28.3    2021-12-09 [2] Bioconductor
 # BiocSingular             1.10.0    2021-10-26 [2] Bioconductor
 # BiocVersion              3.14.0    2021-05-19 [2] Bioconductor
 # biomaRt                  2.50.3    2022-02-03 [2] Bioconductor
 # Biostrings               2.62.0    2021-10-26 [2] Bioconductor
 # bit                      4.0.4     2020-08-04 [2] CRAN (R 4.1.0)
 # bit64                    4.0.5     2020-08-30 [2] CRAN (R 4.1.0)
 # bitops                   1.0-7     2021-04-24 [2] CRAN (R 4.1.0)
 # blob                     1.2.3     2022-04-10 [2] CRAN (R 4.1.2)
 # bluster                  1.4.0     2021-10-26 [2] Bioconductor
 # cachem                   1.0.6     2021-08-19 [2] CRAN (R 4.1.2)
 # callr                    3.7.0     2021-04-20 [2] CRAN (R 4.1.0)
 # cli                      3.3.0     2022-04-25 [2] CRAN (R 4.1.2)
 # cluster                  2.1.3     2022-03-28 [3] CRAN (R 4.1.2)
 # colorspace               2.0-3     2022-02-21 [2] CRAN (R 4.1.2)
 # crayon                   1.5.1     2022-03-26 [2] CRAN (R 4.1.2)
 # curl                     4.3.2     2021-06-23 [2] CRAN (R 4.1.0)
 # DBI                      1.1.3     2022-06-18 [2] CRAN (R 4.1.2)
 # dbplyr                   2.2.1     2022-06-27 [2] CRAN (R 4.1.2)
 # DelayedArray             0.20.0    2021-10-26 [2] Bioconductor
 # DelayedMatrixStats       1.16.0    2021-10-26 [2] Bioconductor
 # devtools               * 2.4.3     2021-11-30 [2] CRAN (R 4.1.2)
 # digest                   0.6.29    2021-12-01 [2] CRAN (R 4.1.2)
 # dplyr                    1.0.9     2022-04-28 [2] CRAN (R 4.1.2)
 # dqrng                    0.3.0     2021-05-01 [2] CRAN (R 4.1.2)
 # DropletUtils           * 1.14.2    2022-01-09 [2] Bioconductor
 # edgeR                    3.36.0    2021-10-26 [2] Bioconductor
 # ellipsis                 0.3.2     2021-04-29 [2] CRAN (R 4.1.0)
 # ensembldb                2.18.4    2022-03-24 [2] Bioconductor
 # ExperimentHub            2.2.1     2022-01-23 [2] Bioconductor
 # fansi                    1.0.3     2022-03-24 [2] CRAN (R 4.1.2)
 # fastmap                  1.1.0     2021-01-25 [2] CRAN (R 4.1.0)
 # filelock                 1.0.2     2018-10-05 [2] CRAN (R 4.1.0)
 # fs                       1.5.2     2021-12-08 [2] CRAN (R 4.1.2)
 # gargle                   1.2.0     2021-07-02 [2] CRAN (R 4.1.0)
 # generics                 0.1.3     2022-07-05 [2] CRAN (R 4.1.2)
 # GenomeInfoDb           * 1.30.1    2022-01-30 [2] Bioconductor
 # GenomeInfoDbData         1.2.7     2021-11-01 [2] Bioconductor
 # GenomicAlignments        1.30.0    2021-10-26 [2] Bioconductor
 # GenomicFeatures          1.46.5    2022-02-27 [2] Bioconductor
 # GenomicRanges          * 1.46.1    2021-11-18 [2] Bioconductor
 # ggbeeswarm               0.6.0     2017-08-07 [2] CRAN (R 4.1.2)
 # ggplot2                * 3.3.6     2022-05-03 [2] CRAN (R 4.1.2)
 # ggrepel                  0.9.1     2021-01-15 [2] CRAN (R 4.1.0)
 # glue                     1.6.2     2022-02-24 [2] CRAN (R 4.1.2)
 # googledrive              2.0.0     2021-07-08 [2] CRAN (R 4.1.0)
 # gridExtra                2.3       2017-09-09 [2] CRAN (R 4.1.0)
 # gtable                   0.3.0     2019-03-25 [2] CRAN (R 4.1.0)
 # HDF5Array                1.22.1    2021-11-14 [2] Bioconductor
 # here                   * 1.0.1     2020-12-13 [2] CRAN (R 4.1.2)
 # hms                      1.1.1     2021-09-26 [2] CRAN (R 4.1.2)
 # htmltools                0.5.2     2021-08-25 [2] CRAN (R 4.1.2)
 # httpuv                   1.6.5     2022-01-05 [2] CRAN (R 4.1.2)
 # httr                     1.4.3     2022-05-04 [2] CRAN (R 4.1.2)
 # igraph                   1.3.2     2022-06-13 [2] CRAN (R 4.1.2)
 # interactiveDisplayBase   1.32.0    2021-10-26 [2] Bioconductor
 # IRanges                * 2.28.0    2021-10-26 [2] Bioconductor
 # irlba                    2.3.5     2021-12-06 [2] CRAN (R 4.1.2)
 # jaffelab               * 0.99.32   2022-03-15 [1] Github (LieberInstitute/jaffelab@7b7afe3)
 # KEGGREST                 1.34.0    2021-10-26 [2] Bioconductor
 # later                    1.3.0     2021-08-18 [2] CRAN (R 4.1.2)
 # lattice                  0.20-45   2021-09-22 [3] CRAN (R 4.1.2)
 # lazyeval                 0.2.2     2019-03-15 [2] CRAN (R 4.1.0)
 # lifecycle                1.0.1     2021-09-24 [2] CRAN (R 4.1.2)
 # limma                    3.50.3    2022-04-07 [2] Bioconductor
 # locfit                   1.5-9.5   2022-03-03 [2] CRAN (R 4.1.2)
 # magrittr                 2.0.3     2022-03-30 [2] CRAN (R 4.1.2)
 # MASS                     7.3-56    2022-03-23 [3] CRAN (R 4.1.2)
 # Matrix                 * 1.4-1     2022-03-23 [3] CRAN (R 4.1.2)
 # MatrixGenerics         * 1.6.0     2021-10-26 [2] Bioconductor
 # matrixStats            * 0.62.0    2022-04-19 [2] CRAN (R 4.1.2)
 # memoise                  2.0.1     2021-11-26 [2] CRAN (R 4.1.2)
 # metapod                  1.2.0     2021-10-26 [2] Bioconductor
 # mime                     0.12      2021-09-28 [2] CRAN (R 4.1.2)
 # munsell                  0.5.0     2018-06-12 [2] CRAN (R 4.1.0)
 # pillar                   1.7.0     2022-02-01 [2] CRAN (R 4.1.2)
 # pkgbuild                 1.3.1     2021-12-20 [2] CRAN (R 4.1.2)
 # pkgconfig                2.0.3     2019-09-22 [2] CRAN (R 4.1.0)
 # pkgload                  1.3.0     2022-06-27 [2] CRAN (R 4.1.2)
 # png                      0.1-7     2013-12-03 [2] CRAN (R 4.1.0)
 # prettyunits              1.1.1     2020-01-24 [2] CRAN (R 4.1.0)
 # processx                 3.7.0     2022-07-07 [2] CRAN (R 4.1.2)
 # progress                 1.2.2     2019-05-16 [2] CRAN (R 4.1.0)
 # promises                 1.2.0.1   2021-02-11 [2] CRAN (R 4.1.0)
 # ProtGenerics             1.26.0    2021-10-26 [2] Bioconductor
 # ps                       1.7.1     2022-06-18 [2] CRAN (R 4.1.2)
 # purrr                    0.3.4     2020-04-17 [2] CRAN (R 4.1.0)
 # R.methodsS3              1.8.2     2022-06-13 [2] CRAN (R 4.1.2)
 # R.oo                     1.25.0    2022-06-12 [2] CRAN (R 4.1.2)
 # R.utils                  2.12.0    2022-06-28 [2] CRAN (R 4.1.2)
 # R6                       2.5.1     2021-08-19 [2] CRAN (R 4.1.2)
 # rafalib                * 1.0.0     2015-08-09 [1] CRAN (R 4.1.2)
 # rappdirs                 0.3.3     2021-01-31 [2] CRAN (R 4.1.0)
 # RColorBrewer             1.1-3     2022-04-03 [2] CRAN (R 4.1.2)
 # Rcpp                     1.0.8.3   2022-03-17 [2] CRAN (R 4.1.2)
 # RCurl                    1.98-1.7  2022-06-09 [2] CRAN (R 4.1.2)
 # remotes                  2.4.2     2021-11-30 [1] CRAN (R 4.1.2)
 # ResidualMatrix           1.4.0     2021-10-26 [1] Bioconductor
 # restfulr                 0.0.15    2022-06-16 [2] CRAN (R 4.1.2)
 # rhdf5                    2.38.1    2022-03-10 [2] Bioconductor
 # rhdf5filters             1.6.0     2021-10-26 [2] Bioconductor
 # Rhdf5lib                 1.16.0    2021-10-26 [2] Bioconductor
 # rjson                    0.2.21    2022-01-09 [2] CRAN (R 4.1.2)
 # rlang                    1.0.3     2022-06-27 [2] CRAN (R 4.1.2)
 # rprojroot                2.0.3     2022-04-02 [2] CRAN (R 4.1.2)
 # Rsamtools                2.10.0    2021-10-26 [2] Bioconductor
 # RSQLite                  2.2.14    2022-05-07 [2] CRAN (R 4.1.2)
 # rsvd                     1.0.5     2021-04-16 [2] CRAN (R 4.1.2)
 # rtracklayer              1.54.0    2021-10-26 [2] Bioconductor
 # Rtsne                  * 0.16      2022-04-17 [2] CRAN (R 4.1.2)
 # S4Vectors              * 0.32.4    2022-03-24 [2] Bioconductor
 # ScaledMatrix             1.2.0     2021-10-26 [2] Bioconductor
 # scales                   1.2.0     2022-04-13 [2] CRAN (R 4.1.2)
 # scater                 * 1.22.0    2021-10-26 [2] Bioconductor
 # scran                  * 1.22.1    2021-11-14 [2] Bioconductor
 # scRNAseq               * 2.8.0     2021-10-30 [1] Bioconductor
 # scuttle                * 1.4.0     2021-10-26 [2] Bioconductor
 # segmented                1.4-1     2022-03-24 [1] CRAN (R 4.1.2)
 # sessioninfo              1.2.2     2021-12-06 [2] CRAN (R 4.1.2)
 # shiny                    1.7.1     2021-10-02 [2] CRAN (R 4.1.2)
 # SingleCellExperiment   * 1.16.0    2021-10-26 [2] Bioconductor
 # sparseMatrixStats        1.6.0     2021-10-26 [2] Bioconductor
 # statmod                  1.4.36    2021-05-10 [2] CRAN (R 4.1.0)
 # stringi                  1.7.6     2021-11-29 [2] CRAN (R 4.1.2)
 # stringr                  1.4.0     2019-02-10 [2] CRAN (R 4.1.0)
 # SummarizedExperiment   * 1.24.0    2021-10-26 [2] Bioconductor
 # tibble                   3.1.7     2022-05-03 [2] CRAN (R 4.1.2)
 # tidyselect               1.1.2     2022-02-21 [2] CRAN (R 4.1.2)
 # usethis                * 2.1.6     2022-05-25 [2] CRAN (R 4.1.2)
 # utf8                     1.2.2     2021-07-24 [2] CRAN (R 4.1.0)
 # uwot                   * 0.1.11    2021-12-02 [2] CRAN (R 4.1.2)
 # vctrs                    0.4.1     2022-04-13 [2] CRAN (R 4.1.2)
 # vipor                    0.4.5     2017-03-22 [2] CRAN (R 4.1.2)
 # viridis                  0.6.2     2021-10-13 [2] CRAN (R 4.1.2)
 # viridisLite              0.4.0     2021-04-13 [2] CRAN (R 4.1.0)
 # withr                    2.5.0     2022-03-03 [2] CRAN (R 4.1.2)
 # XML                      3.99-0.10 2022-06-09 [2] CRAN (R 4.1.2)
 # xml2                     1.3.3     2021-11-30 [2] CRAN (R 4.1.2)
 # xtable                   1.8-4     2019-04-21 [2] CRAN (R 4.1.0)
 # XVector                  0.34.0    2021-10-26 [2] Bioconductor
 # yaml                     2.3.5     2022-02-21 [2] CRAN (R 4.1.2)
 # zlibbioc                 1.40.0    2021-10-26 [2] Bioconductor
 #
 # [1] /users/jstolz/R/4.1.x
 # [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/site-library
 # [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/library
