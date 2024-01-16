#!/bin/bash
#$ -cwd
#$ -l leek,mem_free=3G,h_vmem=5G,h_fsize=150G
#$ -N step00-merge-Roche_Habenula.PairedEnd
#$ -pe local 8
#$ -o ./logs/merge-Roche_Habenula.txt
#$ -e ./logs/merge-Roche_Habenula.txt
#$ -hold_jid pipeline_setup
#$ -m a

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

Rscript /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/step00-merge.R -s /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/.samples_unmerged.manifest -o /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/merged_fastq -c 8

echo "**** Job ends ****"
date
