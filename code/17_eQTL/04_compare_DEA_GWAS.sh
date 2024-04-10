#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=25G
#SBATCH --job-name=04_compare_DEA_GWAS
#SBATCH -o /dev/null
#SBATCH -e /dev/null

run_mode=independent

log_path="logs/04_compare_DEA_GWAS_${run_mode}.log"

{
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
Rscript 04_compare_DEA_GWAS.R -m $run_mode

echo "**** Job ends ****"
date
} > $log_path 2>&1
