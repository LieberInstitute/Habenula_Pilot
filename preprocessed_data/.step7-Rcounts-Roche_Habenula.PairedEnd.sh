#!/bin/bash
#$ -cwd
#$ -pe local 5
#$ -l bluejay,mem_free=18G,h_vmem=20G,h_fsize=200G
#$ -N step7-Rcounts-Roche_Habenula.PairedEnd
#$ -o ./logs/Rcounts-Roche_Habenula.txt
#$ -e ./logs/Rcounts-Roche_Habenula.txt
#$ -hold_jid pipeline_setup,step4-featCounts-Roche_Habenula.PairedEnd,step6-txQuant-Roche_Habenula.PairedEnd
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

## Force R 3.3.x in JHPCE (to avoid some issues with conda_R)
module unload conda_R
module load R/3.3.x

Rscript /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/.step7-create_count_objects-human.R -o hg38 -m /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data -e Roche_Habenula -p PairedEnd -l TRUE -c TRUE -t 5 -s reverse

echo "**** Job ends ****"
date
