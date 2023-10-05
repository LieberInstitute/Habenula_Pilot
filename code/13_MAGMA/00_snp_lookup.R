library("biomaRt")
library("data.table")
library("here")


snp_mart <- useEnsembl(biomart="ENSEMBL_MART_SNP", 
                       host="https://grch37.ensembl.org", 
                       dataset="hsapiens_snp")

## good data
scz_snploc <- fread(here("processed-data", "13_MAGMA","GWAS", "SCZ", "PGC3_SCZ_wave3.european.autosome.public.v3.snploc"))
## data with missing snploc
gwas_mdd <- fread(here("processed-data", "13_MAGMA","GWAS", "mdd2019edinburgh", "PGC_UKB_depression_genome-wide.txt"))
gwas_oud <- fread(here("processed-data", "13_MAGMA","GWAS", "sud2020op", "opi.DEPvEXP_EUR.noAF.tbl"))

 ## unknown 
snps_unknown <- union(gwas_mdd$MarkerName, gwas_oud$rsID)
length(snps_unknown)
# [1] 8876095

snps_known <- scz_snploc[scz_snploc$SNP %in% snps_unknown, ]
dim(snps_known)
# [1] 7333623       3

snps_unknown <- snps_unknown[!snps_unknown %in% snps_known$SNP]
length(snps_unknown)
# [1] 1542472

message(Sys.time(), " - Start snp lookup with bioMart")
snploc <- getBM(attributes = c('refsnp_id', 'chr_name', 'chrom_start'), 
      values = snps_unknown, 
      filters ="snp_filter",
      mart = snp_mart)

message(Sys.time(), " - Saving")
save(snploc, file = here("processed-data", "13_MAGMA","GWAS","snp_mart_lookup.Rdata"))

# slurmjobs::job_single(name = "00_snp_lookup", memory = "25G", create_shell = TRUE, command = "Rscript 00_snp_lookup.R")

load(here("processed-data", "13_MAGMA","GWAS","snp_mart_lookup.Rdata"))
dim(snploc)
# [1] 1600390       3

table(snploc$chr_name)
table(grepl("H",snploc$chr_name))
# FALSE    TRUE 
# 1516456   83934 

snploc <- snploc[!grepl("H",snploc$chr_name),]
dim(snploc)
# [1] 1516456       3

colnames(snploc) <- c("SNP", "CHR", "BP")

snploc <- rbind(snploc, snps_known)
dim(snploc)
rownames(snploc) <- snploc$SNP
# [1] 8850079       3

## MDD
snploc_mdd <- snploc[gwas_mdd$MarkerName,]
table(is.na(snploc_mdd$SNP))
# FALSE    TRUE 
# 8481694    1607
snploc_mdd <- snploc_mdd[!is.na(snploc_mdd$SNP),]

write.table(snploc_mdd,
            file = here("processed-data", "13_MAGMA", "GWAS", "mdd2019edinburgh", "PGC_UKB_depression_genome-wide.snploc"),
            sep = "\t", col.names = T, row.names = F, quote = F
)

## MDD
snploc_oud <- snploc[gwas_oud$rsID,]
table(is.na(snploc_oud$SNP))
# FALSE    TRUE 
# 8481694    1607
snploc_oud <- snploc_oud[!is.na(snploc_oud$SNP),]

write.table(snploc_oud,
            file = here("processed-data", "13_MAGMA", "GWAS", "sud2020op", "opi.DEPvEXP_EUR.noAF.snploc"),
            sep = "\t", col.names = T, row.names = F, quote = F
)
