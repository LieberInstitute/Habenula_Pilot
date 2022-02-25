library("SingleCellExperiment")
library("scran")
library("scater")
library("sessioninfo")
library("here")

load(here("processed-data","08_snRNA-seq_Erik", "20210525_human_hb_processing.rda"), verbose = TRUE)


stats.hb <- perCellQCMetrics(s3e.hb, subsets=list(Mito=grep("^MT-", rowData(s3e.hb)$Symbol)))
high.mito.hb <- isOutlier(stats.hb$subsets_Mito_percent, nmads=3, type="higher")
colData(s3e.hb) <- cbind(colData(s3e.hb), stats.hb, "high.mito" = high.mito.hb)


# Let's also look at low library size and low detected features

#No outliers...so let's set a manual threshold
qc.lib<-s3e.hb$sum < 1000
qc.detected<-s3e.hb$detected<500

s3e.hb$discard <- qc.lib | qc.detected | high.mito.hb

s3e.hb$qc.lib.manual <- qc.lib
s3e.hb$qc.detected.manual <- qc.detected
s3e.hb$sample_short <- basename(gsub("/outs/raw_feature_bc_matrix", "", s3e.hb$Sample))
s3e.hb$qc.lib <- isOutlier(s3e.hb$sum, log=TRUE, type="lower", batch = s3e.hb$Sample, subset = s3e.hb$sample_short != "Br1204")
s3e.hb$qc.detected <- isOutlier(s3e.hb$detected, log=TRUE, type="lower", batch = s3e.hb$Sample, subset = s3e.hb$sample_short != "Br1204")
s3e.hb$high.mito.sample <- isOutlier(stats.hb$subsets_Mito_percent, nmads=3, type="higher", batch = s3e.hb$Sample)


pdf(here("plots","09_snRNA-seq_re-processed","live_checks.pdf"), width = 21)
plotColData(s3e.hb, x = "sample_short", y="sum", colour_by="qc.lib") +
    scale_y_log10() + ggtitle("Total count") + geom_hline(yintercept = 1000)
plotColData(s3e.hb, x = "sample_short", y="detected", colour_by="qc.detected") +
    scale_y_log10() + ggtitle("Detected features") + geom_hline(yintercept = 500)
plotColData(s3e.hb, x = "sample_short", y="subsets_Mito_percent",
              colour_by="high.mito.sample") + geom_hline(yintercept = attr(high.mito.hb, "thresholds")["higher"])
dev.off()


s3e.hb$discard_auto <- s3e.hb$qc.lib | s3e.hb$qc.detected | s3e.hb$high.mito.sample
sce <- s3e.hb

sce$key <- paste0(sce$sample_short, "_", sce$Barcode)

load(here("processed-data","08_snRNA-seq_Erik", "s3e_hb.rda"), verbose = TRUE)
s3e.hb$key <- paste0(s3e.hb$sample_name, "_", s3e.hb$Barcode)

m <- match(s3e.hb$key, sce$key)
s3e.hb$qc.lib <- sce$qc.lib[m]
s3e.hb$qc.detected <- sce$qc.detected[m]
s3e.hb$high.mito.sample <- sce$high.mito.sample[m]
s3e.hb$discard_auto <- sce$discard_auto[m]

addmargins(table("Josh QC" = sce$discard, "Auto QC" = sce$discard_auto))
#        Auto QC
# Josh QC FALSE  TRUE   Sum
#   FALSE 16833   675 17508
#   TRUE    840  3700  4540
#   Sum   17673  4375 22048
addmargins(table(s3e.hb$cellType, s3e.hb$qc.lib))

  #              FALSE  TRUE   Sum
  # Astro          603     6   609
  # Endo           111     7   118
  # Micro          366    13   379
  # Oligo.1       1769    31  1800
  # Oligo.2        426     0   426
  # OPC.1          252    26   278
  # OPC.2          760     5   765
  # Neuron.Ambig    46    33    79
  # LHb.1         1358     0  1358
  # LHb.2          765     1   766
  # LHb.3          615     3   618
  # LHb.4          510     6   516
  # LHb.5          302     0   302
  # LHb.6           73     1    74
  # MHb.1          495     1   496
  # MHb.2          197     0   197
  # Thal.GABA.1   2894    14  2908
  # Thal.GABA.2   4686    15  4701
  # Thal.GABA.3     86     8    94
  # Thal.MD        260     1   261
  # Thal.PF        238     1   239
  # Thal.PVT       545     0   545
  # Sum          17357   172 17529
addmargins(table(s3e.hb$cellType, s3e.hb$qc.detected))
  #              FALSE  TRUE   Sum
  # Astro          592    17   609
  # Endo           100    18   118
  # Micro          359    20   379
  # Oligo.1       1755    45  1800
  # Oligo.2        426     0   426
  # OPC.1          244    34   278
  # OPC.2          756     9   765
  # Neuron.Ambig    42    37    79
  # LHb.1         1358     0  1358
  # LHb.2          764     2   766
  # LHb.3          614     4   618
  # LHb.4          507     9   516
  # LHb.5          302     0   302
  # LHb.6           73     1    74
  # MHb.1          493     3   496
  # MHb.2          197     0   197
  # Thal.GABA.1   2887    21  2908
  # Thal.GABA.2   4679    22  4701
  # Thal.GABA.3     82    12    94
  # Thal.MD        257     4   261
  # Thal.PF        233     6   239
  # Thal.PVT       545     0   545
  # Sum          17265   264 17529
addmargins(table(s3e.hb$cellType, s3e.hb$high.mito.sample))
  #              FALSE  TRUE   Sum
  # Astro          526    83   609
  # Endo            98    20   118
  # Micro          346    33   379
  # Oligo.1       1729    71  1800
  # Oligo.2        425     1   426
  # OPC.1          265    13   278
  # OPC.2          752    13   765
  # Neuron.Ambig    39    40    79
  # LHb.1         1347    11  1358
  # LHb.2          762     4   766
  # LHb.3          618     0   618
  # LHb.4          512     4   516
  # LHb.5          302     0   302
  # LHb.6           74     0    74
  # MHb.1          495     1   496
  # MHb.2          197     0   197
  # Thal.GABA.1   2826    82  2908
  # Thal.GABA.2   4625    76  4701
  # Thal.GABA.3     90     4    94
  # Thal.MD        237    24   261
  # Thal.PF        226    13   239
  # Thal.PVT       543     2   545
  # Sum          17034   495 17529
x <- addmargins(table(s3e.hb$cellType, s3e.hb$discard_auto))
round(100 * sweep(x, 1, x[, 3], "/"), 2)
  #               FALSE   TRUE    Sum
  # Astro         83.74  16.26 100.00
  # Endo          72.88  27.12 100.00
  # Micro         86.54  13.46 100.00
  # Oligo.1       94.17   5.83 100.00
  # Oligo.2       99.77   0.23 100.00
  # OPC.1         84.89  15.11 100.00
  # OPC.2         97.25   2.75 100.00
  # Neuron.Ambig  16.46  83.54 100.00
  # LHb.1         99.19   0.81 100.00
  # LHb.2         99.22   0.78 100.00
  # LHb.3         99.35   0.65 100.00
  # LHb.4         97.48   2.52 100.00
  # LHb.5        100.00   0.00 100.00
  # LHb.6         98.65   1.35 100.00
  # MHb.1         99.19   0.81 100.00
  # MHb.2        100.00   0.00 100.00
  # Thal.GABA.1   96.49   3.51 100.00
  # Thal.GABA.2   98.00   2.00 100.00
  # Thal.GABA.3   84.04  15.96 100.00
  # Thal.MD       90.42   9.58 100.00
  # Thal.PF       92.89   7.11 100.00
  # Thal.PVT      99.63   0.37 100.00
  # Sum           95.94   4.06 100.00

sessioninfo::session_info()

# ─ Session info  ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  hash: sad but relieved face, classical building, woman police officer: medium-light skin tone
#
#  setting  value
#  version  R version 4.1.0 Patched (2021-05-18 r80330)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2022-02-25
#  pandoc   2.12 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1/bin/pandoc
#
# ─ Packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version  date (UTC) lib source
#  assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
#  beachmat               2.8.1    2021-08-12 [2] Bioconductor
#  beeswarm               0.4.0    2021-06-01 [1] CRAN (R 4.1.0)
#  Biobase              * 2.52.0   2021-05-19 [2] Bioconductor
#  BiocGenerics         * 0.38.0   2021-05-19 [2] Bioconductor
#  BiocNeighbors          1.10.0   2021-05-19 [1] Bioconductor
#  BiocParallel           1.26.2   2021-08-22 [2] Bioconductor
#  BiocSingular           1.8.1    2021-06-08 [1] Bioconductor
#  bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.1.0)
#  bluster                1.2.1    2021-05-27 [1] Bioconductor
#  cli                    3.1.0    2021-10-27 [2] CRAN (R 4.1.0)
#  cluster                2.1.2    2021-04-17 [3] CRAN (R 4.1.0)
#  colorspace             2.0-2    2021-06-24 [2] CRAN (R 4.1.0)
#  crayon                 1.4.2    2021-10-29 [2] CRAN (R 4.1.0)
#  DBI                    1.1.1    2021-01-15 [2] CRAN (R 4.1.0)
#  DelayedArray           0.18.0   2021-05-19 [2] Bioconductor
#  DelayedMatrixStats     1.14.3   2021-08-26 [2] Bioconductor
#  digest                 0.6.28   2021-09-23 [2] CRAN (R 4.1.0)
#  dplyr                  1.0.7    2021-06-18 [2] CRAN (R 4.1.0)
#  dqrng                  0.3.0    2021-05-01 [1] CRAN (R 4.1.0)
#  edgeR                  3.34.1   2021-09-05 [2] Bioconductor
#  ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
#  fansi                  0.5.0    2021-05-25 [2] CRAN (R 4.1.0)
#  farver                 2.1.0    2021-02-28 [2] CRAN (R 4.1.0)
#  generics               0.1.1    2021-10-25 [2] CRAN (R 4.1.0)
#  GenomeInfoDb         * 1.28.4   2021-09-05 [2] Bioconductor
#  GenomeInfoDbData       1.2.6    2021-05-11 [2] Bioconductor
#  GenomicRanges        * 1.44.0   2021-05-19 [2] Bioconductor
#  ggbeeswarm             0.6.0    2017-08-07 [1] CRAN (R 4.1.0)
#  ggplot2              * 3.3.5    2021-06-25 [2] CRAN (R 4.1.0)
#  glue                   1.4.2    2020-08-27 [2] CRAN (R 4.1.0)
#  gridExtra              2.3      2017-09-09 [2] CRAN (R 4.1.0)
#  gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
#  here                 * 1.0.1    2020-12-13 [1] CRAN (R 4.1.0)
#  igraph                 1.2.7    2021-10-15 [2] CRAN (R 4.1.0)
#  IRanges              * 2.26.0   2021-05-19 [2] Bioconductor
#  irlba                  2.3.3    2019-02-05 [2] CRAN (R 4.1.0)
#  labeling               0.4.2    2020-10-20 [2] CRAN (R 4.1.0)
#  lattice                0.20-45  2021-09-22 [3] CRAN (R 4.1.0)
#  lifecycle              1.0.1    2021-09-24 [1] CRAN (R 4.1.0)
#  limma                  3.48.3   2021-08-10 [2] Bioconductor
#  locfit                 1.5-9.4  2020-03-25 [2] CRAN (R 4.1.0)
#  magrittr               2.0.1    2020-11-17 [2] CRAN (R 4.1.0)
#  Matrix                 1.3-4    2021-06-01 [3] CRAN (R 4.1.0)
#  MatrixGenerics       * 1.4.3    2021-08-26 [2] Bioconductor
#  matrixStats          * 0.61.0   2021-09-17 [2] CRAN (R 4.1.0)
#  metapod                1.0.0    2021-05-19 [1] Bioconductor
#  munsell                0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
#  pillar                 1.6.3    2021-09-26 [1] CRAN (R 4.1.0)
#  pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
#  purrr                  0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
#  R6                     2.5.1    2021-08-19 [1] CRAN (R 4.1.0)
#  Rcpp                   1.0.7    2021-07-07 [2] CRAN (R 4.1.0)
#  RCurl                  1.98-1.5 2021-09-17 [2] CRAN (R 4.1.0)
#  rlang                  0.4.12   2021-10-18 [2] CRAN (R 4.1.0)
#  rprojroot              2.0.2    2020-11-15 [2] CRAN (R 4.1.0)
#  rsvd                   1.0.5    2021-04-16 [1] CRAN (R 4.1.0)
#  S4Vectors            * 0.30.2   2021-10-03 [2] Bioconductor
#  ScaledMatrix           1.0.0    2021-05-19 [1] Bioconductor
#  scales                 1.1.1    2020-05-11 [2] CRAN (R 4.1.0)
#  scater               * 1.20.1   2021-06-15 [1] Bioconductor
#  scran                * 1.20.1   2021-05-24 [1] Bioconductor
#  scuttle              * 1.2.1    2021-08-05 [1] Bioconductor
#  sessioninfo          * 1.2.1    2021-11-02 [2] CRAN (R 4.1.0)
#  SingleCellExperiment * 1.14.1   2021-05-21 [2] Bioconductor
#  sparseMatrixStats      1.4.2    2021-08-08 [2] Bioconductor
#  statmod                1.4.36   2021-05-10 [2] CRAN (R 4.1.0)
#  SummarizedExperiment * 1.22.0   2021-05-19 [2] Bioconductor
#  tibble                 3.1.5    2021-09-30 [1] CRAN (R 4.1.0)
#  tidyselect             1.1.1    2021-04-30 [2] CRAN (R 4.1.0)
#  utf8                   1.2.2    2021-07-24 [2] CRAN (R 4.1.0)
#  vctrs                  0.3.8    2021-04-29 [2] CRAN (R 4.1.0)
#  vipor                  0.4.5    2017-03-22 [1] CRAN (R 4.1.0)
#  viridis                0.6.2    2021-10-13 [2] CRAN (R 4.1.0)
#  viridisLite            0.4.0    2021-04-13 [2] CRAN (R 4.1.0)
#  withr                  2.4.2    2021-04-18 [2] CRAN (R 4.1.0)
#  XVector                0.32.0   2021-05-19 [2] Bioconductor
#  zlibbioc               1.38.0   2021-05-19 [2] Bioconductor
#
#  [1] /users/jstolz/R/4.1
#  [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1/R/4.1/lib64/R/site-library
#  [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1/R/4.1/lib64/R/library
#
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
