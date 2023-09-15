#!/bin/bash -l
#SBATCH --output=logs/03_MAGMA_gene-set_SCZ.txt
#SBATCH --error=logs/03_MAGMA_gene-set_SCZ.txt
#SBATCH --partition=shared
#SBATCH --job-name=MAGMA_gene-set_SCZ
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

magma --gene-results ../../processed-data/13_MAGMA/SCZ/PGC3_SCZ_wave3.european.autosome.public.v3.ensembl.genes.raw\
	--set-annot ../../processed-data/13_MAGMA/gene_sets/markerSets_broad_ENSEMBL_FDR05.txt gene-col=Gene set-col=Set\
	--out ../../processed-data/13_MAGMA/SCZ_broad

echo "**** Job ends ****"
date
