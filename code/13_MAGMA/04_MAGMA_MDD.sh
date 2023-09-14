#!/bin/bash -l
#SBATCH --output=logs/04_MAGMA_MDD.txt
#SBATCH --error=logs/02_MAGMA_MDD.txt
#SBATCH --partition=shared
#SBATCH --job-name=MAGMA_MDD
#SBATCH --mem=25GB

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## load modules
module load magma/1.10
# list modules
module list

ANNOT_PREFIX="../../processed-data/13_MAGMA/GWAS/MDD/MDD.phs001672.pha005122"

## Step 1 Annotate w/ HGCh38
magma --annotate --snp-loc ../../processed-data/13_MAGMA/GWAS/MDD/MDD.phs001672.pha005122.snploc\
	 --gene-loc ../../processed-data/13_MAGMA/geneloc/GRCh38-ensembl93_to_hg19-lifted_30k-expressing-GENES.gene.loc\
		 --out $ANNOT_PREFIX 

magma --bfile /dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/Analysis/Layer_Guesses/MAGMA/g1000_eur\
	--pval ../../processed-data/13_MAGMA/MDD/MDD.phs001672.pha005122.pval ncol=N\
	--gene-annot $ANNOT_PREFIX.genes.annot\
	--out $ANNOT_PREFIX

magma --gene-results $ANNOT_PREFIX.genes.raw\
	--set-annot ../../processed-data/13_MAGMA/gene_sets/markerSets_broad_ENSEMBL_FDR05.txt gene-col=Gene set-col=Set\
	--out ../../processed-data/13_MAGMA/MAGMA_output/MDD_broad/MDD_braod

echo "**** Job ends ****"
date
