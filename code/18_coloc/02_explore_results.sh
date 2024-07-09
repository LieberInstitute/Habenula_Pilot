#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH --job-name=02_explore_results
#SBATCH -o logs/02_explore_results.log
#SBATCH -e logs/02_explore_results.log

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load conda_R/4.3.x
Rscript 02_explore_results.R

echo "**** Job ends ****"
date
