#!/bin/bash
#$ -cwd
#$ -N step4-featCounts-Roche_Habenula.PairedEnd_clean
#$ -o ./logs/featCounts-Roche_Habenula_clean.txt
#$ -e ./logs/featCounts-Roche_Habenula_clean.txt
#$ -hold_jid pipeline_setup,step4-featCounts-Roche_Habenula.PairedEnd
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

## Delete temporary files after they have been used
rm -rf /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/Counts/junction/tmpdir

echo "**** Job ends ****"
date
