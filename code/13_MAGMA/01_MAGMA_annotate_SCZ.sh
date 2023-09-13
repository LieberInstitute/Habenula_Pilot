#!/bin/bash -l
#SBATCH --output=logs/01_MAGMA_annotate_SCZ.txt
#SBATCH --error=logs/01_MAGMA_annotate_SCZ.txt
#SBATCH --partition=shared
#SBATCH --job-name=MAGMA_annotate_SCZ
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

magma --annotate\
	--snp-loc ../../processed-data/13_MAGMA/SCZ/PGC3_SCZ_wave3.european.autosome.public.v3.snploc\
	--gene-loc ../../processed-data/13_MAGMA/NCBI37.3/NCBI37.3.gene.loc\
	--out ../../processed-data/13_MAGMA/SCZ/PGC3_SCZ_wave3.european.autosome.public.v3

echo "**** Job ends ****"
date