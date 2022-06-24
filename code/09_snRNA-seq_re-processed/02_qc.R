library("SingleCellExperiment")
library("scran")
library("scater")
library("sessioninfo")
library("here")
library(scDblFinder)
library(purrr)

load(here("processed-data","08_snRNA-seq_Erik", "20220301_human_hb_processing.rda"), verbose = TRUE)


stats.hb <- perCellQCMetrics(sce.all.hb, subsets=list(Mito=grep("^MT-", rowData(sce.all.hb)$gene_name)))
high.mito.hb <- isOutlier(stats.hb$subsets_Mito_percent, nmads=3, type="higher")
colData(sce.all.hb ) <- cbind(colData(sce.all.hb), stats.hb, "high.mito" = high.mito.hb)


# Let's also look at low library size and low detected features

#No outliers...so let's set a manual threshold
qc.lib<-sce.all.hb$sum < 1000
qc.detected<-sce.all.hb$detected<500

sce.all.hb$discard <- qc.lib | qc.detected | high.mito.hb

sce.all.hb$qc.lib.manual <- qc.lib
sce.all.hb$qc.detected.manual <- qc.detected
sce.all.hb$qc.lib <- isOutlier(sce.all.hb$sum, log=TRUE, type="lower", batch = sce.all.hb$Sample, subset = sce.all.hb$sample_short != "Br1204")
sce.all.hb$qc.detected <- isOutlier(sce.all.hb$detected, log=TRUE, type="lower", batch = sce.all.hb$Sample, subset = sce.all.hb$sample_short != "Br1204")
sce.all.hb$high.mito.sample <- isOutlier(stats.hb$subsets_Mito_percent, nmads=3, type="higher", batch = sce.all.hb$Sample)


pdf(here("plots","09_snRNA-seq_re-processed","live_checks.pdf"), width = 21)
plotColData(sce.all.hb, x = "sample_short", y="sum", colour_by="qc.lib") +
    scale_y_log10() + ggtitle("Total count") + geom_hline(yintercept = 1000)
plotColData(sce.all.hb, x = "sample_short", y="detected", colour_by="qc.detected") +
    scale_y_log10() + ggtitle("Detected features") + geom_hline(yintercept = 500)
plotColData(sce.all.hb, x = "sample_short", y="subsets_Mito_percent",
              colour_by="high.mito.sample") + geom_hline(yintercept = attr(high.mito.hb, "thresholds")["higher"])
dev.off()


sce.all.hb$discard_auto <- sce.all.hb$qc.lib | sce.all.hb$qc.detected | sce.all.hb$high.mito.sample



load(here("processed-data","08_snRNA-seq_Erik", "s3e_hb.rda"), verbose = TRUE)
s3e.hb$key <- paste0(s3e.hb$sample_name, "_", s3e.hb$Barcode)

m <- match(s3e.hb$key, sce.all.hb$key)
s3e.hb$qc.lib <- sce.all.hb$qc.lib[m]
s3e.hb$qc.detected <- sce.all.hb$qc.detected[m]
s3e.hb$high.mito.sample <- sce.all.hb$high.mito.sample[m]
s3e.hb$discard_auto <- sce.all.hb$discard_auto[m]

addmargins(table("Josh QC" = sce.all.hb$discard, "Auto QC" = sce.all.hb$discard_auto))
#        Auto QC
# Josh QC FALSE  TRUE   Sum
#   FALSE 16066   703 16769
#   TRUE   1050  1983  3033
#   Sum   17116  2686 19802
addmargins(table(sce.all.hb$cellType_Erik, sce.all.hb$qc.lib))
  #
  #              FALSE  TRUE   Sum
  # Astro          602     7   609
  # Endo           110     8   118
  # Micro          366    13   379
  # Oligo.1       1767    32  1799
  # Oligo.2        426     0   426
  # OPC.1          243    29   272
  # OPC.2          759     6   765
  # Neuron.Ambig    45    26    71
  # LHb.1         1357     0  1357
  # LHb.2          749     2   751
  # LHb.3          614     2   616
  # LHb.4          500     5   505
  # LHb.5          302     0   302
  # LHb.6           72     2    74
  # MHb.1          488     3   491
  # MHb.2          197     0   197
  # Thal.GABA.1   2718    13  2731
  # Thal.GABA.2   4344    15  4359
  # Thal.GABA.3     86     8    94
  # Thal.MD        259     1   260
  # Thal.PF        235     2   237
  # Thal.PVT       545     0   545
  # Sum          16784   174 16958
addmargins(table(sce.all.hb$cellType_Erik, sce.all.hb$qc.detected))
  #              FALSE  TRUE   Sum
  # Astro          581    28   609
  # Endo            97    21   118
  # Micro          348    31   379
  # Oligo.1       1740    59  1799
  # Oligo.2        426     0   426
  # OPC.1          236    36   272
  # OPC.2          755    10   765
  # Neuron.Ambig    41    30    71
  # LHb.1         1357     0  1357
  # LHb.2          749     2   751
  # LHb.3          614     2   616
  # LHb.4          500     5   505
  # LHb.5          302     0   302
  # LHb.6           72     2    74
  # MHb.1          487     4   491
  # MHb.2          197     0   197
  # Thal.GABA.1   2713    18  2731
  # Thal.GABA.2   4334    25  4359
  # Thal.GABA.3     82    12    94
  # Thal.MD        256     4   260
  # Thal.PF        231     6   237
  # Thal.PVT       545     0   545
  # Sum          16663   295 16958
addmargins(table(sce.all.hb$cellType_Erik, sce.all.hb$high.mito.sample))

               FALSE  TRUE   Sum
  Astro          581    28   609
  Endo            97    21   118
  Micro          348    31   379
  Oligo.1       1740    59  1799
  Oligo.2        426     0   426
  OPC.1          236    36   272
  OPC.2          755    10   765
  Neuron.Ambig    41    30    71
  LHb.1         1357     0  1357
  LHb.2          749     2   751
  LHb.3          614     2   616
  LHb.4          500     5   505
  LHb.5          302     0   302
  LHb.6           72     2    74
  MHb.1          487     4   491
  MHb.2          197     0   197
  Thal.GABA.1   2713    18  2731
  Thal.GABA.2   4334    25  4359
  Thal.GABA.3     82    12    94
  Thal.MD        256     4   260
  Thal.PF        231     6   237
  Thal.PVT       545     0   545
  Sum          16663   295 16958
> addmargins(table(sce.all.hb$cellType_Erik, sce.all.hb$high.mito.sample))
  #
  #              FALSE  TRUE   Sum
  # Astro          521    88   609
  # Endo            96    22   118
  # Micro          346    33   379
  # Oligo.1       1723    76  1799
  # Oligo.2        425     1   426
  # OPC.1          256    16   272
  # OPC.2          752    13   765
  # Neuron.Ambig    22    49    71
  # LHb.1         1341    16  1357
  # LHb.2          745     6   751
  # LHb.3          615     1   616
  # LHb.4          499     6   505
  # LHb.5          302     0   302
  # LHb.6           74     0    74
  # MHb.1          490     1   491
  # MHb.2          197     0   197
  # Thal.GABA.1   2643    88  2731
  # Thal.GABA.2   4280    79  4359
  # Thal.GABA.3     90     4    94
  # Thal.MD        232    28   260
  # Thal.PF        224    13   237
  # Thal.PVT       541     4   545
  # Sum          16414   544 16958
x <- addmargins(table(sce.all.hb$cellType_Erik, sce.all.hb$discard_auto))
round(100 * sweep(x, 1, x[, 3], "/"), 2)

  #               FALSE   TRUE    Sum
  # Astro         82.10  17.90 100.00
  # Endo          69.49  30.51 100.00
  # Micro         84.43  15.57 100.00
  # Oligo.1       93.61   6.39 100.00
  # Oligo.2       99.77   0.23 100.00
  # OPC.1         82.72  17.28 100.00
  # OPC.2         97.12   2.88 100.00
  # Neuron.Ambig   9.86  90.14 100.00
  # LHb.1         98.82   1.18 100.00
  # LHb.2         98.93   1.07 100.00
  # LHb.3         99.51   0.49 100.00
  # LHb.4         97.82   2.18 100.00
  # LHb.5        100.00   0.00 100.00
  # LHb.6         97.30   2.70 100.00
  # MHb.1         98.98   1.02 100.00
  # MHb.2        100.00   0.00 100.00
  # Thal.GABA.1   96.16   3.84 100.00
  # Thal.GABA.2   97.75   2.25 100.00
  # Thal.GABA.3   84.04  15.96 100.00
  # Thal.MD       88.85  11.15 100.00
  # Thal.PF       92.83   7.17 100.00
  # Thal.PVT      99.27   0.73 100.00
  # Sum           95.48   4.52 100.00

set.seed(506)

colData(sce.all.hb)$doubletScore <- NA

for (i in splitit(sce.all.hb$sample_short)) {
    sce_temp <- sce.all.hb[, i]
    message(ncol(sce_temp))
    ## To speed up, run on sample-level top-HVGs - just take top 1000
    normd <- logNormCounts(sce_temp)
    geneVar <- modelGeneVar(normd)
    topHVGs <- getTopHVGs(geneVar, n = 1000)

    dbl_dens <- computeDoubletDensity(normd, subset.row = topHVGs)
    colData(sce.all.hb)$doubletScore[i] <- dbl_dens
}

save(sce.all.hb, file = here("processed-data","08_snRNA-seq_Erik", "01_qc.rda"))

#  sessioninfo::session_info()
#  package              * version  date (UTC) lib source
#  assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
#  beachmat               2.10.0   2021-10-26 [2] Bioconductor
#  beeswarm               0.4.0    2021-06-01 [2] CRAN (R 4.1.2)
#  Biobase              * 2.54.0   2021-10-26 [2] Bioconductor
#  BiocGenerics         * 0.40.0   2021-10-26 [2] Bioconductor
#  BiocNeighbors          1.12.0   2021-10-26 [2] Bioconductor
#  BiocParallel           1.28.3   2021-12-09 [2] Bioconductor
#  BiocSingular           1.10.0   2021-10-26 [2] Bioconductor
#  bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.1.0)
#  bluster                1.4.0    2021-10-26 [2] Bioconductor
#  cli                    3.3.0    2022-04-25 [2] CRAN (R 4.1.2)
#  cluster                2.1.3    2022-03-28 [3] CRAN (R 4.1.2)
#  colorspace             2.0-3    2022-02-21 [2] CRAN (R 4.1.2)
#  cowplot                1.1.1    2020-12-30 [2] CRAN (R 4.1.2)
#  crayon                 1.5.1    2022-03-26 [2] CRAN (R 4.1.2)
#  data.table             1.14.2   2021-09-27 [2] CRAN (R 4.1.2)
#  DBI                    1.1.2    2021-12-20 [2] CRAN (R 4.1.2)
#  DelayedArray           0.20.0   2021-10-26 [2] Bioconductor
#  DelayedMatrixStats     1.16.0   2021-10-26 [2] Bioconductor
#  digest                 0.6.29   2021-12-01 [2] CRAN (R 4.1.2)
#  dplyr                  1.0.9    2022-04-28 [2] CRAN (R 4.1.2)
#  dqrng                  0.3.0    2021-05-01 [2] CRAN (R 4.1.2)
#  edgeR                  3.36.0   2021-10-26 [2] Bioconductor
#  ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
#  fansi                  1.0.3    2022-03-24 [2] CRAN (R 4.1.2)
#  farver                 2.1.0    2021-02-28 [2] CRAN (R 4.1.0)
#  generics               0.1.2    2022-01-31 [2] CRAN (R 4.1.2)
#  GenomeInfoDb         * 1.30.1   2022-01-30 [2] Bioconductor
#  GenomeInfoDbData       1.2.7    2021-11-01 [2] Bioconductor
#  GenomicRanges        * 1.46.1   2021-11-18 [2] Bioconductor
#  ggbeeswarm             0.6.0    2017-08-07 [2] CRAN (R 4.1.2)
#  ggplot2              * 3.3.6    2022-05-03 [2] CRAN (R 4.1.2)
#  ggrepel                0.9.1    2021-01-15 [2] CRAN (R 4.1.0)
#  glue                   1.6.2    2022-02-24 [2] CRAN (R 4.1.2)
#  gridExtra              2.3      2017-09-09 [2] CRAN (R 4.1.0)
#  gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
#  here                 * 1.0.1    2020-12-13 [2] CRAN (R 4.1.2)
#  igraph                 1.3.2    2022-06-13 [2] CRAN (R 4.1.2)
#  IRanges              * 2.28.0   2021-10-26 [2] Bioconductor
#  irlba                  2.3.5    2021-12-06 [2] CRAN (R 4.1.2)
#  jsonlite               1.8.0    2022-02-22 [2] CRAN (R 4.1.2)
#  labeling               0.4.2    2020-10-20 [2] CRAN (R 4.1.0)
#  lattice                0.20-45  2021-09-22 [3] CRAN (R 4.1.2)
#  lifecycle              1.0.1    2021-09-24 [2] CRAN (R 4.1.2)
#  limma                  3.50.3   2022-04-07 [2] Bioconductor
#  locfit                 1.5-9.5  2022-03-03 [2] CRAN (R 4.1.2)
#  magrittr               2.0.3    2022-03-30 [2] CRAN (R 4.1.2)
#  MASS                   7.3-56   2022-03-23 [3] CRAN (R 4.1.2)
#  Matrix                 1.4-1    2022-03-23 [3] CRAN (R 4.1.2)
#  MatrixGenerics       * 1.6.0    2021-10-26 [2] Bioconductor
#  matrixStats          * 0.62.0   2022-04-19 [2] CRAN (R 4.1.2)
#  metapod                1.2.0    2021-10-26 [2] Bioconductor
#  munsell                0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
#  pillar                 1.7.0    2022-02-01 [2] CRAN (R 4.1.2)
#  pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
#  purrr                * 0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
#  R6                     2.5.1    2021-08-19 [2] CRAN (R 4.1.2)
#  Rcpp                   1.0.8.3  2022-03-17 [2] CRAN (R 4.1.2)
#  RCurl                  1.98-1.7 2022-06-09 [2] CRAN (R 4.1.2)
#  rlang                  1.0.2    2022-03-04 [2] CRAN (R 4.1.2)
#  rprojroot              2.0.3    2022-04-02 [2] CRAN (R 4.1.2)
#  rsvd                   1.0.5    2021-04-16 [2] CRAN (R 4.1.2)
#  S4Vectors            * 0.32.4   2022-03-24 [2] Bioconductor
#  ScaledMatrix           1.2.0    2021-10-26 [2] Bioconductor
#  scales                 1.2.0    2022-04-13 [2] CRAN (R 4.1.2)
#  scater               * 1.22.0   2021-10-26 [2] Bioconductor
#  scDblFinder          * 1.8.0    2021-10-26 [1] Bioconductor
#  scran                * 1.22.1   2021-11-14 [2] Bioconductor
#  scuttle              * 1.4.0    2021-10-26 [2] Bioconductor
#  sessioninfo          * 1.2.2    2021-12-06 [2] CRAN (R 4.1.2)
#  SingleCellExperiment * 1.16.0   2021-10-26 [2] Bioconductor
#  sparseMatrixStats      1.6.0    2021-10-26 [2] Bioconductor
#  statmod                1.4.36   2021-05-10 [2] CRAN (R 4.1.0)
#  SummarizedExperiment * 1.24.0   2021-10-26 [2] Bioconductor
#  tibble                 3.1.7    2022-05-03 [2] CRAN (R 4.1.2)
#  tidyselect             1.1.2    2022-02-21 [2] CRAN (R 4.1.2)
#  utf8                   1.2.2    2021-07-24 [2] CRAN (R 4.1.0)
#  vctrs                  0.4.1    2022-04-13 [2] CRAN (R 4.1.2)
#  vipor                  0.4.5    2017-03-22 [2] CRAN (R 4.1.2)
#  viridis                0.6.2    2021-10-13 [2] CRAN (R 4.1.2)
#  viridisLite            0.4.0    2021-04-13 [2] CRAN (R 4.1.0)
#  withr                  2.5.0    2022-03-03 [2] CRAN (R 4.1.2)
#  xgboost                1.5.2.1  2022-02-21 [1] CRAN (R 4.1.2)
#  XVector                0.34.0   2021-10-26 [2] Bioconductor
#  zlibbioc               1.40.0   2021-10-26 [2] Bioconductor
#
#  [1] /users/jstolz/R/4.1.x
#  [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/site-library
#  [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
