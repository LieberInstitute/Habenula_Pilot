library("biomaRt")
library("data.table")
library("here")


snp_mart <- useEnsembl(biomart="ENSEMBL_MART_SNP", 
                       host="https://grch37.ensembl.org", 
                       dataset="hsapiens_snp")

## good data
scz_snploc <- fread(here("processed-data", "13_MAGMA","GWAS", "SCZ", "PGC3_SCZ_wave3.european.autosome.public.v3.snploc"))

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
      values = snps_unknown[1:100], 
      filters ="snp_filter",
      mart = snp_mart)

save(snploc, file = here(here("processed-data", "13_MAGMA","GWAS","snp_mart_lookup.Rdata")))

# slurmjobs::job_single(name = "00_snp_lookup", memory = "25G", create_shell = TRUE, command = "Rscript 00_snp_lookup.R")

