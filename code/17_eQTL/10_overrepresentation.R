#   The following code was intended to determine if genes or SNPs in significant
#   eQTLs were overrepresented in external reference sets
library(sessioninfo)

#   Columns: gene is vs. is not habenula DEG (FDR < 0.1)
#   Rows: gene is vs. is not in an eQTL (FDR < 0.05)
message("Fisher test for overrepresentation of eQTL genes in habenula DEGs:")
m = matrix(c(7, 166, 710, 21873), ncol = 2)
fisher.test(m, alternative = "greater")

#   Columns: SNP is vs. is not a PGC3 GWAS risk SNP
#   Rows: SNP is vs. is not in an eQTL (FDR < 0.05)
message("Fisher test for overrepresentation of eQTL SNPs in PGC3 GWAS risk SNPs:")
m = matrix(c(15, 18041, 614, 4629710), ncol = 2)
fisher.test(m, alternative = "greater")

session_info()
