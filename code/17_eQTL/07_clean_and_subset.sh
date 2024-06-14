#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=40G
#SBATCH --job-name=07_clean_and_subset
#SBATCH -o logs/07_clean_and_subset.log
#SBATCH -e logs/07_clean_and_subset.log

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load conda_R/4.3.x
Rscript 07_clean_and_subset.R

echo "**** Job ends ****"
date
