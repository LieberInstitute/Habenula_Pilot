#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=50G,h_vmem=50G,h_fsize=100G
#$ -N get_droplet_scores
#$ -o logs/get_droplet_scores.$TASK_ID.txt
#$ -e logs/get_droplet_scores.$TASK_ID.txt
#$ -m e
#$ -t 1-7
#$ -tc 10

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R

## List current modules for reproducibility
module list

## Edit with your job command
samp=$(awk '{ print $0}' /dcl02/lieber/ajaffe/Roche_Habenula/code/09_snRNA-seq_re-processed/droplet_scores_troubleshoot/sceSamples.txt | awk "NR==${SGE_TASK_ID}")

echo $samp

Rscript 01_get_droplet_scores.R $samp 300

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
