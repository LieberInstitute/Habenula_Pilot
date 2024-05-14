#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=20G
#SBATCH --job-name=01_vcf2snpPCs
#SBATCH -c 1
#SBATCH -t 2-00:00:00
#SBATCH -o logs/01_vcf2snpPCs_rerun.txt
#SBATCH -e logs/01_vcf2snpPCs_rerun.txt

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load conda_R/4.3.x
module load plink/1.90b

## List current modules for reproducibility
module list

Rscript 01_vcf2snpPCs.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.1
## available from http://research.libd.org/slurmjobs/
