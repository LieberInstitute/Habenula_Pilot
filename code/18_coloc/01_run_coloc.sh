#!/bin/bash

#SBATCH -p shared
#SBATCH -c 8
#SBATCH --mem=15G
#SBATCH --job-name=01_run_coloc
#SBATCH -o logs/01_run_coloc.log
#SBATCH -e logs/01_run_coloc.log

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
Rscript 01_run_coloc.R

echo "**** Job ends ****"
date
