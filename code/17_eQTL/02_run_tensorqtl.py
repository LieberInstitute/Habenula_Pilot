import my_tensorqtl_run
from tensorqtl import cis
import sys
import pandas as pd
import session_info

import pandas as pd
import torch
import tensorqtl
from tensorqtl import genotypeio, cis, trans
print(f'PyTorch {torch.__version__}')
print(f'Pandas {pd.__version__}')

def load_data(plink_prefix_path, expression_bed, covariates_file, add_chr = False, fix_geno_names = False):
    
    # load phenotypes and covariates
    print("Reading Expression files: " + expression_bed)
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
    
    ## Fix genoSample names
    if fix_geno_names:
        print("Using fam to assign genoSample names")
        fam = pd.read_table(plink_prefix_path + '.fam', delimiter = " ", names = ["V" + str(i) for i in range(6)])
        genoSamples = [str(v1) + "_" + str(v0) for (v0,v1) in zip(fam.V1, fam.V0)]
        genotype_df.columns = genoSamples
    
    ## Fix chr names (maybe fix in plink?)
    if add_chr:
        print("Adding 'chr' to genotype positions")
        variant_df.chrom = [s.split(':')[0] for s in list(variant_df.index)]

    # Filter expression to chrom in snp data
    variant_chrom = set(variant_df.chrom)
    express_chrom = set(phenotype_pos_df.chr)
    
    if express_chrom - variant_chrom:
        print("Excluding phenotypes from these chromosomes:")
        print(express_chrom - variant_chrom)
    
        chrom_filter = phenotype_pos_df.chr.isin(variant_chrom)
        print("Phenotypes with chr in variants: " )
        print(chrom_filter.value_counts())
    
        phenotype_df = phenotype_df[chrom_filter]
        phenotype_pos_df = phenotype_pos_df[chrom_filter]
    
    return(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df)
    
def get_input_paths(feature, region):
    expression_bed = '../data/tensorQTL_input/expression_bed/' + feature + '_' + region + '.bed.gz'
    covariates_file = '../data/tensorQTL_input/covariates_txt/covariates_' + feature + '_' + region + '.txt'
    return expression_bed, covariates_file

pair = sys.argv[1].split("_")
feature = pair[0]
region = pair[1]

print("FEATURE = " + feature)
print("REGION = " + region)

expres, covar = my_tensorqtl_run.get_input_paths(feature, region)

print("Expression bed file: " + expres)
print("Covar file: " + covar)

prefix = '/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/eqtl/data/tensorQTL_out/genomewide_nominal'
plink = "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/genotype_data/mdd_bpd/maf01/mdd_bpd_maf01"

genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df = my_tensorqtl_run.load_data(plink, expres, covar, add_chr = True, fix_geno_names = True)

tag = prefix +"/" + feature + "_" + region
print('\n**** STARTING tensorQTL ****')
print("Saving output to: " + tag)

cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, prefix = tag, covariates_df= covariates_df,
                maf_threshold=0.05, interaction_df=None, maf_threshold_interaction=0, 
                group_s=None, window=500000, run_eigenmt=True, output_dir= tag, write_top=False, verbose=False)

session_info.show()