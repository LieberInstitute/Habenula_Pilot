library("SingleCellExperiment")
library("scran")
library("scater")
library("sessioninfo")
library("here")

load(here("processed-data","08_snRNA-seq_Erik", "20220301_human_hb_processing.rda"), verbose = TRUE)


stats.hb <- perCellQCMetrics(sce.all.hb, subsets=list(Mito=grep("^MT-", rowData(sce.all.hb)$Symbol)))
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
#   FALSE 17897   329 18226
#   TRUE    391  1247  1638
#   Sum   18288  1576 19864
addmargins(table(sce.all.hb$cellType_Erik, sce.all.hb$qc.lib))
  #
  #              FALSE  TRUE   Sum
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
addmargins(table(sce.all.hb$cellType_Erik, sce.all.hb$qc.detected))
  #               FALSE  TRUE   Sum
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
addmargins(table(sce.all.hb$cellType_Erik, sce.all.hb$high.mito.sample))
  #               FALSE   Sum
  # Astro          591   591
  # Endo           110   110
  # Micro          371   371
  # Oligo.1       1783  1783
  # Oligo.2        425   425
  # OPC.1          252   252
  # OPC.2          756   756
  # Neuron.Ambig    46    46
  # LHb.1         1253  1253
  # LHb.2          686   686
  # LHb.3          557   557
  # LHb.4          420   420
  # LHb.5          291   291
  # LHb.6           72    72
  # MHb.1          448   448
  # MHb.2           93    93
  # Thal.GABA.1   2768  2768
  # Thal.GABA.2   4493  4493
  # Thal.GABA.3     87    87
  # Thal.MD        251   251
  # Thal.PF        235   235
  # Thal.PVT       537   537
  # Sum          16525 16525
x <- addmargins(table(sce.all.hb$cellType_Erik, sce.all.hb$discard_auto))
round(100 * sweep(x, 1, x[, 3], "/"), 2)
  #               FALSE   TRUE    Sum
  # Astro         98.82   1.18 100.00
  # Endo          92.73   7.27 100.00
  # Micro         97.04   2.96 100.00
  # Oligo.1       97.87   2.13 100.00
  # Oligo.2      100.00   0.00 100.00
  # OPC.1         88.49  11.51 100.00
  # OPC.2         98.94   1.06 100.00
  # Neuron.Ambig  65.22  34.78 100.00
  # LHb.1         99.84   0.16 100.00
  # LHb.2         99.85   0.15 100.00
  # LHb.3         99.82   0.18 100.00
  # LHb.4         99.05   0.95 100.00
  # LHb.5        100.00   0.00 100.00
  # LHb.6         94.44   5.56 100.00
  # MHb.1         98.88   1.12 100.00
  # MHb.2        100.00   0.00 100.00
  # Thal.GABA.1   99.24   0.76 100.00
  # Thal.GABA.2   99.49   0.51 100.00
  # Thal.GABA.3   87.36  12.64 100.00
  # Thal.MD       98.41   1.59 100.00
  # Thal.PF       97.87   2.13 100.00
  # Thal.PVT     100.00   0.00 100.00
  # Sum           98.80   1.20 100.00

library(scDblFinder)
library(purrr)
sample_id_names<- names(table(sce.all.hb$sample_short))
names(sample_id_names)<-sample_id_names

sample_id_rse<- map(sample_id_names,~sce.all.hb[,sce.all.hb$sample_short==.x])
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
#  0.007  0.050  0.086  0.114  0.143  0.179  0.207  0.236  0.272  0.308  0.351
#    55%    60%    65%    70%    75%    80%    85%    90%    95%   100%
#  0.386  0.429  0.479  0.531  0.622  0.723  0.870  1.095  1.611 31.871
#
# $Br1204
#    0%    5%   10%   15%   20%   25%   30%   35%   40%   45%   50%   55%   60%
# 0.174 0.290 0.353 0.403 0.447 0.507 0.577 0.674 0.770 0.824 0.894 0.983 1.048
#   65%   70%   75%   80%   85%   90%   95%  100%
# 1.127 1.207 1.353 1.470 1.605 1.680 1.848 5.435
#
# $Br1469
#     0%     5%    10%    15%    20%    25%    30%    35%    40%    45%    50%
#  0.005  0.044  0.082  0.115  0.142  0.170  0.197  0.228  0.263  0.304  0.351
#    55%    60%    65%    70%    75%    80%    85%    90%    95%   100%
#  0.397  0.455  0.537  0.625  0.734  0.866  1.086  1.299  1.777 11.448
#
# $Br1735
#     0%     5%    10%    15%    20%    25%    30%    35%    40%    45%    50%
#  0.000  0.019  0.075  0.159  0.215  0.262  0.299  0.346  0.383  0.430  0.476
#    55%    60%    65%    70%    75%    80%    85%    90%    95%   100%
#  0.532  0.598  0.673  0.747  0.859  0.981  1.158  1.420  1.859 14.602
#
# $Br5555
#     0%     5%    10%    15%    20%    25%    30%    35%    40%    45%    50%
#  0.008  0.053  0.084  0.099  0.122  0.137  0.160  0.175  0.198  0.213  0.236
#    55%    60%    65%    70%    75%    80%    85%    90%    95%   100%
#  0.259  0.290  0.320  0.358  0.404  0.473  0.572  0.755  1.090 26.425
#
# $Br5558
#    0%    5%   10%   15%   20%   25%   30%   35%   40%   45%   50%   55%   60%
# 0.044 0.125 0.264 0.340 0.438 0.535 0.629 0.704 0.773 0.859 0.940 1.027 1.127
#   65%   70%   75%   80%   85%   90%   95%  100%
# 1.278 1.449 1.613 1.778 1.956 2.234 2.562 4.804
#
# $Br5639
#     0%     5%    10%    15%    20%    25%    30%    35%    40%    45%    50%
#  0.007  0.048  0.068  0.102  0.129  0.170  0.224  0.272  0.333  0.401  0.483
#    55%    60%    65%    70%    75%    80%    85%    90%    95%   100%
#  0.578  0.680  0.796  0.932  1.082  1.231  1.449  1.762  2.299 17.257
map(dbl.dens.focused,~table(.x >= 5))
# $Br1092
#
# FALSE  TRUE
#  3540    37
#
# $Br1204
#
# FALSE  TRUE
#   713     1
#
# $Br1469
#
# FALSE  TRUE
#  2713    27
#
# $Br1735
#
# FALSE  TRUE
#  4628    43
#
# $Br5555
#
# FALSE  TRUE
#  3737    75
#
# $Br5558
#
# FALSE
#   949
#
# $Br5639
#
# FALSE  TRUE
#  3371    30

# Percent that would be dropped at density score >= 5
round(sapply(names(dbl.dens.focused), function(x) {
  table(dbl.dens.focused[[x]] >= 5)["TRUE"] / ncol(sample_id_rse[[x]]) * 100
}), 3)

# Br1092.TRUE Br1204.TRUE Br1469.TRUE Br1735.TRUE Br5555.TRUE   Br5558.NA
#       1.034       0.140       0.985       0.921       1.967          NA
# Br5639.TRUE
#       0.882

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
