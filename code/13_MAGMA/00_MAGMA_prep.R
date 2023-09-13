
# library("GenomicRanges")
# library("rtracklayer")
library("tidyverse")
library("data.table")
library("here")
library("sessioninfo")

#### GWAS SZC Data ####
gwas_scz = fread(here("processed-data", "13_MAGMA", "SCZ", "PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz"))
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

#### Step 3 Gene-set Analysis ####
GS_PREFIX

"magma --bfile [REFDATA] --pval [PVAL_FILE] N=[N] --gene-annot [ANNOT_PREFIX].genes.annot \ --out [GENE_PREFIX]" 
## create SET_FILE with each row corresponding to a gene set: name of the gene set followed by the gene IDs, separated by whitespace).



# sgejobs::job_single('08_explore_proportions', create_shell = TRUE, queue= 'bluejay', memory = '25G', command = "Rscript 08_explore_proportions.R")



