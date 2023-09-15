# https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/51d15ef9f5f2c4c53f55e22e3fe467de1a724668/MAGMA/10x_setUpMAGMA_MNT2021.R#L12

library(rtracklayer)
library(here)

gtf = import("/dcl01/ajaffe/data/lab/singleCell/refdata-cellranger-hg19-3.0.0/genes/genes.gtf")
## of length 2565061
gtf = gtf[gtf$type == "gene"]
length(gtf)
# [1] 32738

table(seqnames(gtf) %in% c(1:21,"X","Y","MT"))
# FALSE  TRUE 
# 732 32006

gtf <- gtf[seqnames(gtf) %in% c(1:21,"X","Y","MT"),]

geneloc <- data.frame(ensembl = gtf$gene_id, 
                      chr = seqnames(gtf),
                      start = start(gtf),
                      end = end(gtf),
                      strand = strand(gtf),
                      symbol = gtf$gene_name)

head(geneloc)
# ensembl chr  start    end strand       symbol
# 1 ENSG00000243485   1  29554  31109      +   MIR1302-10
# 2 ENSG00000237613   1  34554  36081      -      FAM138A
# 3 ENSG00000186092   1  69091  70008      +        OR4F5
# 4 ENSG00000238009   1  89295 133566      - RP11-34P13.7
# 5 ENSG00000239945   1  89551  91105      - RP11-34P13.8
# 6 ENSG00000237683   1 134901 139379      -   AL627309.1

table(geneloc$chr)
dim(geneloc)
# [1] 32006     6
write.table(annoTab, 
            file= here("processed-data", "13_MAGMA", "geneloc", "GRCh37_Ensembl.gene.loc"), 
            sep="\t",row.names=F, col.names=F, quote=F)
