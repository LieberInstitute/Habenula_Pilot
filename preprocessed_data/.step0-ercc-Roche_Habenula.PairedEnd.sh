#!/bin/bash
#$ -cwd
#$ -l leek,mem_free=3G,h_vmem=5G,h_fsize=100G
#$ -N step0-ercc-Roche_Habenula.PairedEnd
#$ -pe local 8
#$ -o ./logs/ercc-Roche_Habenula.$TASK_ID.txt
#$ -e ./logs/ercc-Roche_Habenula.$TASK_ID.txt
#$ -t 1-73
#$ -tc 20
#$ -hold_jid pipeline_setup,step00-merge-Roche_Habenula.PairedEnd
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"
echo "****"
echo "Sample id: $(cat /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/samples.manifest | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")"
echo "****"

## Locate file and ids
FILE1=$(awk 'BEGIN {FS="\t"} {print $1}' /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/samples.manifest | awk "NR==${SGE_TASK_ID}")
if [ TRUE == "TRUE" ] 
then
    FILE2=$(awk 'BEGIN {FS="\t"} {print $3}' /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/samples.manifest | awk "NR==${SGE_TASK_ID}")
fi
ID=$(cat /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/samples.manifest | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")
mkdir -p /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/Ercc/${ID}

if [ TRUE == "TRUE" ]
then 
    /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/kallisto quant     -i /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/ERCC/ERCC92.idx     -o /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/Ercc/${ID} -t 8 --rf-stranded     ${FILE1} ${FILE2}
else
    /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/kallisto quant     -i /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/ERCC/ERCC92.idx     -o /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/Ercc/${ID} -t 8 --single --rf-stranded ${FILE1}
fi

echo "**** Job ends ****"
date
