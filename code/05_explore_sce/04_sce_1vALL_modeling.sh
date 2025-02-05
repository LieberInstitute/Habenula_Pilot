#!/bin/bash
#$ -cwd
#$ -l mem_free=30G,h_vmem=30G,h_fsize=100G
#$ -N sce_1vALL_modeling
#$ -o logs/04_sce_1vALL_modeling.txt
#$ -e logs/04_sce_1vALL_modeling.txt
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
module load conda_R/4.3

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 04_sce_1vALL_modeling.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
