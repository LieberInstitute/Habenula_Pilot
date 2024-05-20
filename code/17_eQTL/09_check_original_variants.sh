#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=20G
#SBATCH --job-name=09_check_original_variants
#SBATCH -o logs/09_check_original_variants.log
#SBATCH -e logs/09_check_original_variants.log

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load tensorqtl/1.0.8
python 09_check_original_variants.py

echo "**** Job ends ****"
date
