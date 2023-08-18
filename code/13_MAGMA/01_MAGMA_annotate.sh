#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=25G,h_vmem=25G,h_fsize=100G
#$ -N MAGMA_annotate
#$ -o logs/01_MAGMA_annotate.txt
#$ -e logs/01_MAGMA_annotate.txt
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
module load magma/1.10

## List current modules for reproducibility
module list

## Edit with your job command
magma --annotate --snp-loc /dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/08_bulk_snpPC/habenula_genotypes.bim --gene-loc /dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/13_MAGMA/NCBI38/NCBI38.gene.loc --out /dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/13_MAGMA/01_MAGMA_annotate/Habenula_MAGMA
echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
