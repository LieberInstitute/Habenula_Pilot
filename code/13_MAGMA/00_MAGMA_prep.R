
library("tidyverse")
library("data.table")
library("here")
library("sessioninfo")
library("org.Hs.eg.db")


#### GWAS SZC Data ####
gwas = fread(here("processed-data", "13_MAGMA","GWAS", "SCZ", "PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz"))
dim(gwas)
# [1] 7659767      14
head(gwas)
#    CHROM         ID       POS A1 A2  FCAS  FCON IMPINFO         BETA     SE   PVAL  NCAS  NCON     NEFF
# 1:     8 rs62513865 101592213  C  T 0.930 0.927   0.963  0.011997738 0.0171 0.4847 53386 77258 58749.13
# 2:     8 rs79643588 106973048  G  A 0.907 0.906   0.997 -0.008596847 0.0148 0.5605 53386 77258 58749.13
# 3:     8 rs17396518 108690829  T  G 0.565 0.566   0.985 -0.002102208 0.0087 0.8145 53386 77258 58749.13
# 4:     8   rs983166 108681675  A  C 0.564 0.563   0.988  0.004897985 0.0087 0.5704 53386 77258 58749.13
# 5:     8 rs28842593 103044620  T  C 0.840 0.840   0.948 -0.003897586 0.0121 0.7488 53386 77258 58749.13
# 6:     8  rs7014597 104152280  G  C 0.841 0.838   0.994  0.007898723 0.0117 0.5034 53386 77258 58749.13

## snploc
snploc <- gwas |>
  select(SNP = ID, CHR = CHROM, BP = POS)
# only autosomes
snploc |> count(CHR)

write.table(snploc,
            file = here("processed-data", "13_MAGMA", "SCZ", "PGC3_SCZ_wave3.european.autosome.public.v3.snploc"),
            sep = "\t", col.names = T, row.names = F, quote = F
)

table(gwas$NCAS)

## SNP p-vals
snp_pval <- gwas |>
  mutate(N = NCAS + NCON) |>
  select(SNP = ID, P=PVAL, N)

head(snp_pval_scz)

write.table(snp_pval,
            file = here("processed-data", "13_MAGMA", "SCZ", "PGC3_SCZ_wave3.european.autosome.public.v3.pval"),
            sep = "\t", col.names = T, row.names = F, quote = F
)

#### PGC mdd2019edinburgh Data ####
gwas = fread(here("processed-data", "13_MAGMA","GWAS", "mdd2019edinburgh", "PGC_UKB_depression_genome-wide.txt"))
dim(gwas)
# [1] 8483301       7

## no snp location info - snploc created in 00_snp_lookup.R
head(gwas) 
#    MarkerName A1 A2   Freq   LogOR StdErrLogOR         P
# 1:  rs2326918  a  g 0.8452  0.0106      0.0060 0.0756100
# 2:  rs7929618  c  g 0.1314 -0.0224      0.0064 0.0004804
# 3: rs66941928  t  c 0.8031  0.0003      0.0055 0.9502000
# 4:  rs7190157  a  c 0.3517  0.0024      0.0045 0.5992000
# 5: rs12364336  a  g 0.8685  0.0075      0.0064 0.2450000
# 6:  rs6977693  t  c 0.8544  0.0089      0.0061 0.1442000

## SNP p-vals
snp_pval <- gwas |>
  dplyr::select(SNP =  MarkerName, P)

write.table(snp_pval,
            file = here("processed-data", "13_MAGMA", "GWAS", "mdd2019edinburgh", "PGC_UKB_depression_genome-wide.pval"),
            sep = "\t", col.names = T, row.names = F, quote = F
)

#### PGC panic2019 ####
gwas = fread(here("processed-data", "13_MAGMA","GWAS", "panic2019", "pgc-panic2019.vcf.tsv.gz"))
dim(gwas)
# [1] 10151624       16

head(gwas) 
#      #CHROM     POS              ID A1 A2        BETA     SE    PVAL NGT   FCAS   FCON IMPINFO NEFFDIV2 NCAS NCON   DIRE
#   1:     10 1689546      rs11250701  A  G  0.02409731 0.0387 0.53280   0 0.6570 0.6530   0.945  3332.53 2147 7760 +-++++
#   2:     10 2622752 chr10_2622752_D I2  D  0.13370013 0.1191 0.26150   0 0.9760 0.9710   0.944  3332.53 2147 7760 ---+-+
#   3:     10  151476       rs7085086  A  G -0.04210407 0.0398 0.28970   0 0.3080 0.3120   0.949  3332.53 2147 7760 -+--+-
#   4:     10 1593759     rs113494187  T  G  0.33939653 0.1741 0.05117   0 0.9870 0.9840   0.899  3332.53 2147 7760 +++++?
#   5:     10 1708106     rs117915320  A  C -0.39580195 0.2409 0.10030   0 0.0111 0.0132   0.628  3332.53 2147 7760 ?--?--
#   6:     10  790310     rs182753344  T  C -0.02470261 0.0769 0.74840   0 0.0971 0.0950   0.618  3332.53 2147 7760 ---++-

## snploc
snploc <- gwas |>
  dplyr::select(SNP = ID, CHR = `#CHROM`, BP = POS)

head(snploc)

write.table(snploc,
            file = here("processed-data", "13_MAGMA", "GWAS", "panic2019", "pgc-panic2019.snploc"),
            sep = "\t", col.names = T, row.names = F, quote = F
)

## SNP p-vals
snp_pval <- gwas |>
  mutate(N = NCAS + NCON) |>
  dplyr::select(SNP = ID, P=PVAL, N)

#                 SNP       P    N
# 1:       rs11250701 0.53280 9907
# 2:  chr10_2622752_D 0.26150 9907
# 3:        rs7085086 0.28970 9907
# 4:      rs113494187 0.05117 9907
# 5:      rs117915320 0.10030 9907

write.table(snp_pval,
            file = here("processed-data", "13_MAGMA", "GWAS", "panic2019", "pgc-panic2019.pval"),
            sep = "\t", col.names = T, row.names = F, quote = F
)

#### PGC sud2020op ####
# unzipped "OD_cases_vs._opioid-exposed_controls_in_European-ancestry_cohorts.gz"
gwas <- fread(here("processed-data", "13_MAGMA","GWAS", "sud2020op", "opi.DEPvEXP_EUR.noAF.tbl"))
dim(gwas)
# [1] 4211587      14

## no snploc data 
head(gwas) 
#          rsID Allele1 Allele2  Weight Zscore P-value HetISq HetChiSq HetDf HetPVal Total_N Total_NCase Total_NControl
# 1: rs10868284       a       c 3340.40 -1.609  0.1075    0.0    3.923     7 0.78860    5709        3174           2535
# 2: rs62291883       t       c 3483.72 -1.583  0.1134  -10.0    9.998     9 0.35060    6072        3218           2854
# 3: rs61956327       t       c 3043.32 -1.455  0.1458    0.0    4.324     8 0.82680    5260        3081           2179
# 4: rs60994383       a       c 3089.74  1.070  0.2846  -10.9    7.214     6 0.30150    5140        3104           2036
# 5: rs12531896       t       g 3512.19  1.323  0.1857   40.5   20.155    10 0.02782    6115        3262           2853
# 6: rs35515951       a       t 3488.59  1.191  0.2337    0.0    5.705     9 0.76910    6062        3213           2849
#    ngt
# 1:   0
# 2:   0
# 3:   0
# 4:   1
# 5:   1
# 6:   0


## SNP p-vals
snp_pval <- gwas |>
  dplyr::select(SNP = rsID, P=`P-value`, N = Total_N)

#           SNP      P    N
# 1: rs10868284 0.1075 5709
# 2: rs62291883 0.1134 6072
# 3: rs61956327 0.1458 5260
# 4: rs60994383 0.2846 5140
# 5: rs12531896 0.1857 6115

summary(snp_pval$N)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 4137    5571    5876    5759    6050    6148

write.table(snp_pval,
            file = here("processed-data", "13_MAGMA", "GWAS", "sud2020op", "opi.DEPvEXP_EUR.noAF.pval"),
            sep = "\t", col.names = T, row.names = F, quote = F
)


#### GWAS MDD Data ####
gwas_mdd = fread(here("processed-data", "13_MAGMA","GWAS", 'MDD', "MDD.phs001672.pha005122.txt"), skip=21)
dim(gwas_mdd)
# [1] 11700    22
head(gwas_mdd)
# ID Analysis ID     SNP ID   P-value  Rank Plot data Chr ID Chr Position Submitted SNP ID ss2rs rs2genome
# 1: 516388416        5122 rs12137092 4.145e-05  7897         5      1      4666965               ss    NA         +
#   2: 516387459        5122  rs3003457 6.309e-05  9777         5      1     17181435               ss    NA         +
#   3: 516387255        5122   rs754171 8.449e-05 10909         5      1     17196325               ss    NA         +
#   4: 516387447        5122  rs2501818 4.033e-05  7766         5      1     17192406               ss    NA         +
#   5: 516387457        5122  rs2977232 5.965e-05  9552         5      1     17183742               ss    NA         +
#   6: 516387461        5122  rs3104441 5.213e-05  9036         5      1     17191335               ss    NA         -
#   Allele1 Allele2 Minor allele pHWE Call Rate &beta;     SE R-Squared Coded Allele Sample size Bin ID
# 1:       A       G            A   NA        NA 0.0234 0.0057        NA            A          NA      3
# 2:       T       C            T   NA        NA 0.0249 0.0062        NA            T          NA      9
# 3:       T       C            T   NA        NA 0.0244 0.0062        NA            T          NA      9
# 4:       C       G            C   NA        NA 0.0248 0.0060        NA            C          NA      9
# 5:       A       G            G   NA        NA 0.0246 0.0061        NA            A          NA      9
# 6:       A       C            C   NA        NA 0.0245 0.0060        NA            A          NA      9

## snploc
snploc_mdd <- gwas_mdd |>
  dplyr::select(SNP = `SNP ID`, CHR = `Chr ID`, BP = `Chr Position`)

# only autosomes
snploc_mdd |> count(CHR)

write.table(snploc_mdd,
            file = here("processed-data", "13_MAGMA", "GWAS", "MDD", "MDD.phs001672.pha005122.snploc"),
            sep = "\t", col.names = T, row.names = F, quote = F
)

## SNP p-vals
gwas_mdd |> count(`Sample size`)

snp_pval_mdd <- gwas_mdd |>
  dplyr::select(SNP = `SNP ID`, P=`P-value`)

head(snp_pval_mdd)

write.table(snp_pval_mdd,
            file = here("processed-data", "13_MAGMA", 'GWAS',"MDD", "MDD.phs001672.pha005122.pval"),
            sep = "\t", col.names = T, row.names = F, quote = F
)


prep_magma_files <- function(input_file, skip = 21){
  #### GWAS MDD Data ####
  gwas = fread(input_file, skip=21)
  stopifnot(all(c("SNP", "Chr ID", "Chr Position", "P-value") %in% colnames(gwas)))
  ## snploc
  snploc <- gwas |>
    dplyr::select(SNP = `SNP ID`, CHR = `Chr ID`, BP = `Chr Position`)
  
  write.table(snploc,
              file = gsub(".txt$", ".snploc", input_file),
              sep = "\t", col.names = T, row.names = F, quote = F
  )
  
  snp_pval <- gwas |>
    dplyr::select(SNP = `SNP ID`, P=`P-value`)
  
  head(snp_pval)
  
  write.table(snp_pval,
              file = gsub(".txt$", "..pval", input_file),
              sep = "\t", col.names = T, row.names = F, quote = F
  )
}

###OUD data ####
prep_magma_files(here("processed-data", "13_MAGMA", 'GWAS',"OUD", "OUD.phs001672.pha004954.txt"), skip = 17)


#### Step 3 Gene-set Analysis ####

## 1vAll genes from snRNA-seq
load(here("processed-data","05_explore_sce","04_sce_1vALL_modeling","sce_modeling_broad_Annotations.Rdata"), verbose = TRUE)
head(sce_modeling_broad_Annotations$enrichment)

## convert to gene IDs
entrez <- select(org.Hs.eg.db, sce_modeling_broad_Annotations$enrichment$ensembl, 
       columns="ENTREZID", keytype="ENSEMBL") |>
  group_by(ENSEMBL) |>
  dplyr::slice(1) |> # not perfect 1:1...
  dplyr::rename(ensembl = ENSEMBL)

enrichment_long <- sce_modeling_broad_Annotations$enrichment |>
  left_join(entrez) |>
  dplyr::select(gene, ENTREZID, ensembl, starts_with("fdr")) |>
  pivot_longer(!c(gene, ensembl, ENTREZID), names_to = "Set", values_to = 'FDR', names_prefix = "fdr_") |>
  filter(FDR < 0.05) 

enrichment_long |> group_by(Set) |> summarize(n = n(), entrez = sum(is.na(ENTREZID)))
# Set            n entrez
# <chr>      <int>  <int>
# 1 Astrocyte    763    143
# 2 Endo        5125    630
# 3 Excit.Thal   322     98
# 4 Inhib.Thal   512    140
# 5 LHb           25      9
# 6 MHb          127     37
# 7 Microglia   5959    712
# 8 OPC          196     45
# 9 Oligo        269     56


enrichment_long |>
  arrange(Set) |>
  dplyr::select(Set, Gene = ENTREZID) |>
  write.table(file = here("processed-data", "13_MAGMA", "gene_sets", "markerSets_broad_ENTREZID_FDR05.txt"),
              sep = "\t", col.names = T, row.names = F, quote = F
  )

enrichment_long |>
  arrange(Set) |>
  dplyr::select(Set, Gene = ensembl) |>
  write.table(file = here("processed-data", "13_MAGMA", "gene_sets", "markerSets_broad_ENSEMBL_FDR05.txt"),
              sep = "\t", col.names = T, row.names = F, quote = F
  )

