#!/bin/bash -l
#SBATCH --output=logs/03_MAGMA_scz2022.txt
#SBATCH --error=logs/03_MAGMA_scz2022.txt
#SBATCH --partition=shared
#SBATCH --job-name=03_MAGMA_scz2022
#SBATCH --mem=25GB

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

ANNOT_PREFIX="../../processed-data/13_MAGMA/GWAS/scz2022/PGC3_SCZ_wave3.european.autosome.public.v3"

## Step 1 Annotate w/ HGCh38
magma --annotate --snp-loc $ANNOT_PREFIX.snploc\
	--gene-loc ../../processed-data/13_MAGMA/geneloc/GRCh38-ensembl93_to_hg19-lifted_30k-expressing-GENES.gene.loc\
	--out $ANNOT_PREFIX

## Step 2 Gene Analysis
magma --bfile /dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/Analysis/Layer_Guesses/MAGMA/g1000_eur\
	--pval $ANNOT_PREFIX.pval ncol=N\
	--gene-annot $ANNOT_PREFIX.genes.annot\
	--out $ANNOT_PREFIX

## Step 3 Gene Set Analsysis
magma --gene-results $ANNOT_PREFIX.genes.raw\
	--set-annot ../../processed-data/13_MAGMA/gene_sets/markerSets_broad_ENSEMBL_FDR05.txt gene-col=Gene set-col=Set\
	--out ../../processed-data/13_MAGMA/MAGMA_output/scz2022/scz2022_broad

echo "**** Job ends ****"
date
