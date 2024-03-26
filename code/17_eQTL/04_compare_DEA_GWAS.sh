#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=25G
#SBATCH --job-name=04_compare_DEA_GWAS
#SBATCH -o logs/04_compare_DEA_GWAS.log
#SBATCH -e logs/04_compare_DEA_GWAS.log

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load conda_R/4.3.x
module load liftover/1.0
Rscript 04_compare_DEA_GWAS.R

echo "**** Job ends ****"
date
