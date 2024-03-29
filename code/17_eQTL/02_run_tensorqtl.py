import sys
import pandas as pd
import session_info
from pyhere import here
from pathlib import Path

import pandas as pd
import torch
import tensorqtl
from tensorqtl import genotypeio, cis
print(f'PyTorch {torch.__version__}')
print(f'Pandas {pd.__version__}')

run_mode = sys.argv[1]
assert run_mode in ['nominal', 'cis', 'independent', 'interaction']

in_dir = Path(here("processed-data", "17_eQTL", "tensorQTL_input"))
out_dir = Path(here("processed-data", "17_eQTL", "tensorQTL_output", run_mode))
plink_prefix_path = str(
    here("processed-data", '08_bulk_snpPC', "habenula_genotypes")
)
covariates_file = str(in_dir / "covariates.txt")
expression_bed = str(in_dir / "logcounts.bed.gz")

prefix = "habenula"
add_chr = False

out_dir.mkdir(exist_ok = True)

################################################################################
#   Load and ensure compatibility in input data
################################################################################

# load phenotypes and covariates
print("Reading Expression file: " + expression_bed)
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
print("Phenotype dimensions:")
print(phenotype_df.shape)

covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T

print("." *20 )
# PLINK reader for genotypes
print("Reading Plink files: " + plink_prefix_path)
pr = genotypeio.PlinkReader(plink_prefix_path)
print("Loading Genotypes...", end='')
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]
print("Genotype dimensions:", end='')
print(genotype_df.shape)

genotype_ids = set(genotype_df.columns)
covariates_ids = set(covariates_df.index)

#   Exclude any samples that have genotypes but not covariates; halt in the
#   reverse situation
if (genotype_ids > covariates_ids):
    extra_donors = list(genotype_ids - covariates_ids)
    print("Excluding extra genotyped donors:", ','.join(extra_donors))
    genotype_df.drop(extra_donors, axis = 1, inplace = True)
elif (genotype_ids != covariates_ids):
    print("Not all covariate samples in genotyped samples.")
    sys.exit()

## Fix chr names (maybe fix in plink?)
if add_chr:
    print("Adding 'chr' to genotype positions")
    variant_df.chrom = [s.split(':')[0] for s in list(variant_df.index)]

variant_chrom = set(variant_df.chrom)
express_chrom = set(phenotype_pos_df.chr)

#   Exclude any expressed chromosomes that aren't in genotyping data; halt in
#   the reverse situation
if express_chrom > variant_chrom:
    print("Excluding phenotypes from these chromosomes:")
    print(express_chrom - variant_chrom)
    
    chrom_filter = phenotype_pos_df.chr.isin(variant_chrom)
    print("Phenotypes with chr in variants: " )
    print(chrom_filter.value_counts())
    
    phenotype_df = phenotype_df[chrom_filter]
    phenotype_pos_df = phenotype_pos_df[chrom_filter]
elif express_chrom != variant_chrom:
    print("Expression data must contain the same or a superset of chromosomes as the genotyping data.")
    sys.exit()

################################################################################
#   Run tensorQTL
################################################################################

if run_mode == "nominal":
    # nominal code here
elif run_mode == "cis":
    # cis code here
elif run_mode == "independent":
    # independent code here
else:
    #   'run_mode' must be 'interaction' based on check at the top of script

    #   interaction code here

session_info.show()
