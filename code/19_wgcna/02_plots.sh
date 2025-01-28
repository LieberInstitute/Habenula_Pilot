#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=10G
#SBATCH --job-name=02_plots
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o logs/02_plots.txt
#SBATCH -e logs/02_plots.txt

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load the R module
module load conda_R/4.3.x

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 02_plots.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.1
## available from http://research.libd.org/slurmjobs/
