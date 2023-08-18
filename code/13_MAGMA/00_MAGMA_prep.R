
# library("GenomicRanges")
# library("rtracklayer")
# library("tidyverse")
library("here")
library("sessioninfo")


#### SNP loc ####
## use bim file
# "use a .bim file from binary PLINK data, this can be provided to MAGMA as a SNP location file without modification."
# ^ So that was a lie * this bim file didn't work (I think because of 'chr' prefix)

bim_fn <- here("processed-data","08_bulk_snpPC", "habenula_genotypes.bim")
bim <- read.table(bim_fn)
head(bim)
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

snploc <- bim |>
  select(SNP = V2, CHR = V1, BP = V4) |>
  mutate(CHR = gsub("chr","", CHR))

head(snploc)

snploc |> count(CHR)

snploc_fn <- here("processed-data", "13_MAGMA", "habenula_genotypes.snploc")

write.table(snploc,
            file = snploc_fn,
            sep = "\t", col.names = T, row.names = F, quote = F
)

#### Step 1 annotation  ####

## use gene loc build 38 from MAGAM website
geneloc_fn <- here("processed-data", "13_MAGMA", "NCBI38","NCBI38.gene.loc")

output_dir <- here("processed-data", "13_MAGMA", "MAGMA_output")
if(!dir.exists(output_dir)) dir.create(output_dir)

ANNOT_PREFIX <- here(output_dir, "Habenula_MAGMA")

sgejobs::job_single('01_MAGMA_annotate', create_shell = TRUE, memory = '25G', 
                    command = paste0("magma --annotate --snp-loc ", bim_fn,
                                     " --gene-loc ", snploc_fn,
                                     " --out ", ANNOT_PREFIX))

#### Step 2 Gene analysis ####
GENE_PREFIX <- here(output_dir, "Habenula_MAGMA_gene") 

sgejobs::job_single('02_MAGMA_gene', create_shell = TRUE, memory = '25G', 
                    command = paste0("magma --bfile ", "/dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/Analysis/Layer_Guesses/MAGMA/g1000_eur", 
                                     "--gene-annot ", ANNOT_PREFIX,".genes.annot",
                                     "--out ", GENE_PREFIX))


#### Step 3 Gene-set Analysis ####
GS_PREFIX

"magma --bfile [REFDATA] --pval [PVAL_FILE] N=[N] --gene-annot [ANNOT_PREFIX].genes.annot \ --out [GENE_PREFIX]" 
## create SET_FILE with each row corresponding to a gene set: name of the gene set followed by the gene IDs, separated by whitespace).



# sgejobs::job_single('08_explore_proportions', create_shell = TRUE, queue= 'bluejay', memory = '25G', command = "Rscript 08_explore_proportions.R")



