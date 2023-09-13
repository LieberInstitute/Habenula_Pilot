#!/bin/bash -l
#SBATCH --output=logs/02_MAGMA_gene_SCZ.txt
#SBATCH --error=logs/02_MAGMA_gene_SCZ.txt
#SBATCH --partition=shared
#SBATCH --job-name=MAGMA_gene_SCZ
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

magma --bfile /dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/Analysis/Layer_Guesses/MAGMA/g1000_eur\
	--pval ../../processed-data/13_MAGMA/SCZ/PGC3_SCZ_wave3.european.autosome.public.v3.pval ncol=N\
	--gene-annot ../../processed-data/13_MAGMA/SCZ/PGC3_SCZ_wave3.european.autosome.public.v3.genes.annot\
	--out ../../processed-data/13_MAGMA/SCZ/PGC3_SCZ_wave3.european.autosome.public.v3

echo "**** Job ends ****"
date
