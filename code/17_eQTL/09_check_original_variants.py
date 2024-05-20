#   On the first eQTL analysis, 742 independent eQTLs were found to be
#   significant at FDR < 0.05. Later, we fixed the input genotyping data. This
#   script checks if those 742 eQTLs are even in the newly filtered genotyping
#   data, giving an initial idea of how much eQTL results may change due to
#   the genotyping change.

import pandas as pd
import session_info
from pyhere import here
import pandas as pd
from tensorqtl import genotypeio

plink_prefix_path = "/dcs05/lieber/liebercentral/libdGenotype_LIBD001/BrainGenotyping/subsets/Habenula_n69/plink/Hb_gt.hwe_1e-6_geno_05_maf_05"

sig_eqtl_path = here(
    'processed-data', '17_eQTL', 'pre_geno_snapshot_2024_05_02',
    'tensorQTL_output', 'independent', 'FDR05.csv'
)
overlap_eqtl_path = here(
    'processed-data', '17_eQTL', 'pre_geno_snapshot_2024_05_02',
    'rsID_independent_deg_or_gwas_wide.csv'
)

pr = genotypeio.PlinkReader(plink_prefix_path)
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

sig_eqtl = pd.read_csv(sig_eqtl_path)
overlap_eqtl = pd.read_csv(overlap_eqtl_path)

print('Number of originally found significant independent variants measured in newly filtered genotyping data:')
print(sig_eqtl.loc[:, 'variant_id'].isin(variant_df.index).value_counts())

print('Number of originally found variants that overlap GWAS or DEA measured in newly filtered genotyping data:')
print(overlap_eqtl.loc[:, 'variant_id'].isin(variant_df.index).value_counts())

print('Info about the overlapping eQTL that doesnt have a genotyped variant:')
print(overlap_eqtl.loc[~overlap_eqtl.loc[:, 'variant_id'].isin(variant_df.index)])

session_info.show()
