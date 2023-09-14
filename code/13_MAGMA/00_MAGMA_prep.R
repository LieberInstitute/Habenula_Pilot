
library("tidyverse")
library("data.table")
library("here")
library("sessioninfo")
library("org.Hs.eg.db")

#### GWAS SZC Data ####
gwas_scz = fread(here("processed-data", "13_MAGMA","GWAS", "SCZ", "PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz"))
dim(gwas_scz)
# [1] 7659767      14
head(gwas_scz)
#    CHROM         ID       POS A1 A2  FCAS  FCON IMPINFO         BETA     SE   PVAL  NCAS  NCON     NEFF
# 1:     8 rs62513865 101592213  C  T 0.930 0.927   0.963  0.011997738 0.0171 0.4847 53386 77258 58749.13
# 2:     8 rs79643588 106973048  G  A 0.907 0.906   0.997 -0.008596847 0.0148 0.5605 53386 77258 58749.13
# 3:     8 rs17396518 108690829  T  G 0.565 0.566   0.985 -0.002102208 0.0087 0.8145 53386 77258 58749.13
# 4:     8   rs983166 108681675  A  C 0.564 0.563   0.988  0.004897985 0.0087 0.5704 53386 77258 58749.13
# 5:     8 rs28842593 103044620  T  C 0.840 0.840   0.948 -0.003897586 0.0121 0.7488 53386 77258 58749.13
# 6:     8  rs7014597 104152280  G  C 0.841 0.838   0.994  0.007898723 0.0117 0.5034 53386 77258 58749.13

## snploc
snploc_scz <- gwas_scz |>
  select(SNP = ID, CHR = CHROM, BP = POS)
# only autosomes
snploc_scz |> count(CHR)

write.table(snploc_scz,
            file = here("processed-data", "13_MAGMA", "SCZ", "PGC3_SCZ_wave3.european.autosome.public.v3.snploc"),
            sep = "\t", col.names = T, row.names = F, quote = F
)

table(gwas_scz$NCAS)

## SNP p-vals
snp_pval_scz <- gwas_scz |>
  mutate(N = NCAS + NCON) |>
  select(SNP = ID, P=PVAL, N)

head(snp_pval_scz)

write.table(snp_pval_scz,
            file = here("processed-data", "13_MAGMA", "SCZ", "PGC3_SCZ_wave3.european.autosome.public.v3.pval"),
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

