This directory contains genotyping files used as input to the eQTL analysis
performed for this project. The genotyping data was updated twice, so the
`v1` and `v2` directories contain older, unused files. `v3` contains the latest genotyping
data-- in particular, temporary files created while computing
`v3/habenula_R.9_MAF.05.RSann_filt.mds`, the SNP-based PCs. The
`habenula_maf05.bed`, `habenula_maf05.bim`, and `habenula_maf05.fam` Plink
files were used as input to tensorQTL [in this script](https://github.com/LieberInstitute/Habenula_Pilot/blob/master/code/17_eQTL/02_run_tensorqtl.py).
