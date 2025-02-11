#!/bin/bash

#SBATCH -p katun
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH --job-name=10_overrepresentation
#SBATCH -o logs/10_overrepresentation.log
#SBATCH -e logs/10_overrepresentation.log

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load conda_R/4.3.x
Rscript 10_overrepresentation.R

echo "**** Job ends ****"
date
