#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=25G,h_vmem=25G,h_fsize=100G
#$ -N Pseudobulking
#$ -o logs/07c_Pseudobulking.txt
#$ -e logs/07c_Pseudobulking.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R/4.2.x

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 07c_Pseudobulking.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
