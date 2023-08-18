
# library("GenomicRanges")
# library("rtracklayer")
# library("tidyverse")
library("here")
library("sessioninfo")


#### SNP loc ####
## use bim file
# use a .bim file from binary PLINK data, this can be provided to MAGMA as a SNP location file without modification.
bim_fn <- here("processed-data","08_bulk_snpPC", "habenula_genotypes.bim")
# bim <- read.table(bim_fn)
# head(bim)
# # V1                 V2 V3    V4 V5    V6
# # 1 chr1 chr1:11171:CCTTG:C  0 11171  C CCTTG
# # 2 chr1     chr1:23308:G:C  0 23308  C     G
# # 3 chr1    chr1:24963:GT:G  0 24963  G    GT
# # 4 chr1     chr1:48824:T:C  0 48824  C     T
# # 5 chr1     chr1:54490:G:A  0 54490  A     G
# # 6 chr1     chr1:60351:A:G  0 60351  G     A
# 
# nrow(bim)
# [1] 12393872

#### Step 1 annotation  ####

## use gene loc build 38 from MAGAM website
geneloc_fn <- here("processed-data", "13_MAGMA", "NCBI38","NCBI38.gene.loc")

output_dir <- here("processed-data", "13_MAGMA", "01_MAGMA_annotate")

sgejobs::job_single('01_MAGMA_annotate', create_shell = TRUE, memory = '25G', 
                    command = paste0("magma --annotate --snp-loc ", bim_fn,
                                     " --gene-loc ", geneloc_fn,
                                     " --out ", output_dir))

#### Step 2 ####


#### Step 3 Gene-set Analysis ####

## create SET_FILE with each row corresponding to a gene set: name of the gene set followed by the gene IDs, separated by whitespace).



# sgejobs::job_single('08_explore_proportions', create_shell = TRUE, queue= 'bluejay', memory = '25G', command = "Rscript 08_explore_proportions.R")



