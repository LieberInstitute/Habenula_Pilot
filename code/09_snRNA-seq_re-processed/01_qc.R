library("SingleCellExperiment")
library("scran")
library("scater")
library("sessioninfo")
library("here")

load(here("processed-data","08_snRNA-seq_Erik", "20220301_human_hb_processing.rda"), verbose = TRUE)


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
#   FALSE 15891   587 16478
#   TRUE   1005  2381  3386
#   Sum   16896  2968 19864
addmargins(table(s3e.hb$cellType, s3e.hb$qc.lib))
  #           FALSE  TRUE   Sum
  # Astro          585     6   591
  # Endo           102     8   110
  # Micro          360    11   371
  # Oligo.1       1749    34  1783
  # Oligo.2        425     0   425
  # OPC.1          226    26   252
  # OPC.2          750     6   756
  # Neuron.Ambig    32    14    46
  # LHb.1         1251     2  1253
  # LHb.2          685     1   686
  # LHb.3          556     1   557
  # LHb.4          416     4   420
  # LHb.5          291     0   291
  # LHb.6           70     2    72
  # MHb.1          445     3   448
  # MHb.2           93     0    93
  # Thal.GABA.1   2754    14  2768
  # Thal.GABA.2   4476    17  4493
  # Thal.GABA.3     80     7    87
  # Thal.MD        250     1   251
  # Thal.PF        233     2   235
  # Thal.PVT       537     0   537
  # Sum          16366   159 16525
addmargins(table(s3e.hb$cellType, s3e.hb$qc.detected))
  #            FALSE  TRUE   Sum
  # Astro          584     7   591
  # Endo           102     8   110
  # Micro          360    11   371
  # Oligo.1       1745    38  1783
  # Oligo.2        425     0   425
  # OPC.1          223    29   252
  # OPC.2          748     8   756
  # Neuron.Ambig    30    16    46
  # LHb.1         1252     1  1253
  # LHb.2          685     1   686
  # LHb.3          556     1   557
  # LHb.4          416     4   420
  # LHb.5          291     0   291
  # LHb.6           68     4    72
  # MHb.1          443     5   448
  # MHb.2           93     0    93
  # Thal.GABA.1   2747    21  2768
  # Thal.GABA.2   4471    22  4493
  # Thal.GABA.3     76    11    87
  # Thal.MD        247     4   251
  # Thal.PF        230     5   235
  # Thal.PVT       537     0   537
  # Sum          16329   196 16525
addmargins(table(s3e.hb$cellType, s3e.hb$high.mito.sample))
  #           FALSE  TRUE   Sum
  # Astro          521    70   591
  # Endo           100    10   110
  # Micro          343    28   371
  # Oligo.1       1731    52  1783
  # Oligo.2        425     0   425
  # OPC.1          240    12   252
  # OPC.2          755     1   756
  # Neuron.Ambig     7    39    46
  # LHb.1         1235    18  1253
  # LHb.2          678     8   686
  # LHb.3          554     3   557
  # LHb.4          414     6   420
  # LHb.5          291     0   291
  # LHb.6           72     0    72
  # MHb.1          447     1   448
  # MHb.2           93     0    93
  # Thal.GABA.1   2683    85  2768
  # Thal.GABA.2   4427    66  4493
  # Thal.GABA.3     83     4    87
  # Thal.MD        232    19   251
  # Thal.PF        222    13   235
  # Thal.PVT       528     9   537
  # Sum          16081   444 16525
x <- addmargins(table(s3e.hb$cellType, s3e.hb$discard_auto))
round(100 * sweep(x, 1, x[, 3], "/"), 2)

  #               FALSE   TRUE    Sum
  # Astro         86.97  13.03 100.00
  # Endo          88.18  11.82 100.00
  # Micro         89.76  10.24 100.00
  # Oligo.1       95.46   4.54 100.00
  # Oligo.2      100.00   0.00 100.00
  # OPC.1         85.71  14.29 100.00
  # OPC.2         98.94   1.06 100.00
  # Neuron.Ambig   2.17  97.83 100.00
  # LHb.1         98.40   1.60 100.00
  # LHb.2         98.69   1.31 100.00
  # LHb.3         99.28   0.72 100.00
  # LHb.4         97.62   2.38 100.00
  # LHb.5        100.00   0.00 100.00
  # LHb.6         94.44   5.56 100.00
  # MHb.1         98.66   1.34 100.00
  # MHb.2        100.00   0.00 100.00
  # Thal.GABA.1   96.35   3.65 100.00
  # Thal.GABA.2   98.11   1.89 100.00
  # Thal.GABA.3   83.91  16.09 100.00
  # Thal.MD       92.03   7.97 100.00
  # Thal.PF       93.19   6.81 100.00
  # Thal.PVT      98.32   1.68 100.00
  # Sum           96.39   3.61 100.00

library(scDblFinder)
library(purrr)
sample_id_names<- names(table(s3e.hb$sample_name))
names(sample_id_names)<-sample_id_names

sample_id_rse<- map(sample_id_names,~s3e.hb[,s3e.hb$sample_name==.x])
## To speed up, run on sample-level top-HVGs - just take top 1000 ===
pilot.data.normd <- map(sample_id_rse, ~logNormCounts(.x))
geneVar.samples <- map(pilot.data.normd, ~modelGeneVar(.x))
topHVGs <- map(geneVar.samples, ~getTopHVGs(.x, n = 1000))

# Generate doublet density scores
set.seed(109)
dbl.dens.focused <- map(names(pilot.data.normd), ~computeDoubletDensity(pilot.data.normd[[.x]], subset.row=topHVGs[[.x]]))
names(dbl.dens.focused) <- names(pilot.data.normd)

map(dbl.dens.focused, ~round(quantile(.x, probs=seq(0,1,by=0.05)),3))
# $Br1092
#     0%     5%    10%    15%    20%    25%    30%    35%    40%    45%    50%
#  0.007  0.069  0.097  0.125  0.152  0.187  0.228  0.263  0.298  0.332  0.374
#    55%    60%    65%    70%    75%    80%    85%    90%    95%   100%
#  0.415  0.471  0.526  0.602  0.713  0.838  1.018  1.321  1.883 32.640
#
# $Br1204
#    0%    5%   10%   15%   20%   25%   30%   35%   40%   45%   50%   55%   60%
# 0.024 0.120 0.152 0.176 0.205 0.240 0.296 0.341 0.389 0.440 0.493 0.544 0.598
#   65%   70%   75%   80%   85%   90%   95%  100%
# 0.658 0.730 0.810 0.898 1.008 1.183 1.478 6.948
#
# $Br1469
#     0%     5%    10%    15%    20%    25%    30%    35%    40%    45%    50%
#  0.000  0.038  0.066  0.097  0.127  0.156  0.184  0.231  0.283  0.333  0.387
#    55%    60%    65%    70%    75%    80%    85%    90%    95%   100%
#  0.462  0.538  0.627  0.759  0.909  1.122  1.344  1.680  2.206 12.422
#
# $Br1735
#     0%     5%    10%    15%    20%    25%    30%    35%    40%    45%    50%
#  0.000  0.021  0.036  0.057  0.078  0.107  0.149  0.206  0.270  0.348  0.434
#    55%    60%    65%    70%    75%    80%    85%    90%    95%   100%
#  0.540  0.633  0.747  0.889  1.042  1.217  1.450  1.778  2.382 19.752
#
# $Br5555
#     0%     5%    10%    15%    20%    25%    30%    35%    40%    45%    50%
#  0.007  0.052  0.075  0.097  0.112  0.127  0.150  0.165  0.187  0.210  0.240
#    55%    60%    65%    70%    75%    80%    85%    90%    95%   100%
#  0.270  0.300  0.345  0.397  0.472  0.556  0.689  0.876  1.287 28.926
#
# $Br5558
#    0%    5%   10%   15%   20%   25%   30%   35%   40%   45%   50%   55%   60%
# 0.106 0.251 0.370 0.476 0.581 0.664 0.750 0.843 0.934 1.043 1.144 1.277 1.431
#   65%   70%   75%   80%   85%   90%   95%  100%
# 1.646 1.849 2.078 2.279 2.464 2.747 3.406 5.519
#
# $Br5639
#     0%     5%    10%    15%    20%    25%    30%    35%    40%    45%    50%
#  0.000  0.033  0.053  0.084  0.128  0.177  0.226  0.279  0.340  0.416  0.500
#    55%    60%    65%    70%    75%    80%    85%    90%    95%   100%
#  0.579  0.685  0.792  0.902  1.048  1.203  1.406  1.672  2.145 15.840
map(dbl.dens.focused,~table(.x >= 5))
# $Br1092
#
# FALSE  TRUE
#  3426    36
#
# $Br1204
#
# FALSE  TRUE
#  1331     2
#
# $Br1469
#
# FALSE  TRUE
#  2329    29
#
# $Br1735
#
# FALSE  TRUE
#  3509    46
#
# $Br5555
#
# FALSE  TRUE
#  3671    74
#
# $Br5558
#
# FALSE  TRUE
#   862     3
#
# $Br5639
#
# FALSE  TRUE
#  2191    20

# Percent that would be dropped at density score >= 5
round(sapply(names(dbl.dens.focused), function(x) {
  table(dbl.dens.focused[[x]] >= 5)["TRUE"] / ncol(sample_id_rse[[x]]) * 100
}), 3)

# Br1092.TRUE Br1204.TRUE Br1469.TRUE Br1735.TRUE Br5555.TRUE Br5558.TRUE
#       1.040       0.150       1.230       1.294       1.976       0.347
# Br5639.TRUE
#       0.905
# Add the doublet density scores to the colData
for(i in names(sample_id_rse)){
  sample_id_rse[[i]]$doubletScore <- dbl.dens.focused[[i]]
}




 sessioninfo::session_info()
# ─ Session info  ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  hash: men holding hands: medium-dark skin tone, flag: Palestinian Territories, shushing face
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
#  date     2022-03-02
#  pandoc   2.12 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1/bin/pandoc
#
# ─ Packages ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
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
#  data.table             1.14.2   2021-09-27 [2] CRAN (R 4.1.0)
#  DBI                    1.1.1    2021-01-15 [2] CRAN (R 4.1.0)
#  DelayedArray           0.18.0   2021-05-19 [2] Bioconductor
#  DelayedMatrixStats     1.14.3   2021-08-26 [2] Bioconductor
#  dplyr                  1.0.7    2021-06-18 [2] CRAN (R 4.1.0)
#  dqrng                  0.3.0    2021-05-01 [1] CRAN (R 4.1.0)
#  edgeR                  3.34.1   2021-09-05 [2] Bioconductor
#  ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
#  fansi                  0.5.0    2021-05-25 [2] CRAN (R 4.1.0)
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
#  jsonlite               1.7.2    2020-12-09 [2] CRAN (R 4.1.0)
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
#  purrr                * 0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
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
#  scDblFinder          * 1.6.0    2021-05-19 [1] Bioconductor
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
#  vipor                  0.4.5    2017-03-22 [1] C
