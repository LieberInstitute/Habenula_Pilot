#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=3G
#SBATCH --job-name=05_subset_vcf
#SBATCH -o logs/05_subset_vcf.log
#SBATCH -e logs/05_subset_vcf.log

#   The original VCF is too large to read in with VariantAnnotation::readVcf(),
#   but the original method in 04_compare_DEA_GWAS.R (reading in the plink bed)
#   gives different genotype numbers (1-3, which don't obviously correspond to
#   0-2 from readVcf()). This script subsets the big VCF to just the
#   eQTL-significant SNPs whose paired gene is differentially expressed between
#   cases and controls, so VariantAnnotation::readVcf() becomes feasible and can
#   be compared with the current method to verify it works.

big_vcf=../../processed-data/08_bulk_snpPC/habenula_genotypes.vcf.gz
out_dir=../../processed-data/17_eQTL

paired_variants=${out_dir}/DEA_paired_variants.txt
small_vcf=${out_dir}/paired_DEA_eQTL_SNPs.vcf.gz

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

zcat $big_vcf \
    | grep -E "^#|$(paste -sd "|" $paired_variants)" \
    | gzip \
    > $small_vcf
echo "**** Job ends ****"
date
