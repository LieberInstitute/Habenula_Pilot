#!/bin/bash -l
#SBATCH --output=logs/10_MAGMA_gene_set_combo.txt
#SBATCH --error=logs/10_MAGMA_gene_set_combo.txt
#SBATCH --partition=shared
#SBATCH --job-name=10_MAGMA_gene_set_combo
#SBATCH --mem=5GB

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## load modules
module load magma/1.10
# list modules
module list

## mdd2019edinburgh
magma --gene-results ../../processed-data/13_MAGMA/GWAS/mdd2019edinburgh/PGC_UKB_depression_genome-wide\
	--set-annot ../../processed-data/13_MAGMA/gene_sets/markerSets_broad_combo_ENSEMBL_FDR05.txt gene-col=Gene set-col=Set\
	--out ../../processed-data/13_MAGMA/MAGMA_output/mdd2019edinburgh/mdd2019edinburgh_broad_combo

## panic2019
magma --gene-results ../../processed-data/13_MAGMA/GWAS/panic2019/pgc-panic2019\
	--set-annot ../../processed-data/13_MAGMA/gene_sets/markerSets_broad_combo_ENSEMBL_FDR05.txt gene-col=Gene set-col=Set\
	--out ../../processed-data/13_MAGMA/MAGMA_output/panic2019/panic2019_broad_combo

## scz2022
magma --gene-results ../../processed-data/13_MAGMA/GWAS/scz2022/PGC3_SCZ_wave3.european.autosome.public\
	--set-annot ../../processed-data/13_MAGMA/gene_sets/markerSets_broad_combo_ENSEMBL_FDR05.txt gene-col=Gene set-col=Set\
	--out ../../processed-data/13_MAGMA/MAGMA_output/scz2022/scz2022_broad_combo

## sud2020op
magma --gene-results ../../processed-data/13_MAGMA/GWAS/sud2020op/opi.DEPvEXP_EUR.noAF\
	--set-annot ../../processed-data/13_MAGMA/gene_sets/markerSets_broad_combo_ENSEMBL_FDR05.txt gene-col=Gene set-col=Set\
	--out ../../processed-data/13_MAGMA/MAGMA_output/sud2020op/sud2020op_broad_combo


echo "**** Job ends ****"
date
