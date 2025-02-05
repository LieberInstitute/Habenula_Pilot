#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=12G,h_vmem=16G,h_fsize=200G
#$ -N step5-coverage-Roche_Habenula.PairedEnd
#$ -o ./logs/coverage-Roche_Habenula.$TASK_ID.txt
#$ -e ./logs/coverage-Roche_Habenula.$TASK_ID.txt
#$ -t 1-73
#$ -tc 40
#$ -hold_jid pipeline_setup,step3-hisat2-Roche_Habenula.PairedEnd,step3b-infer-strandness-Roche_Habenula.PairedEnd
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

## Force R 3.3.x in JHPCE (to avoid some issues with conda_R)
module unload conda_R
module load R/3.3.x

if [ ! -f "inferred_strandness_pattern.txt" ]
then
    echo "Missing the file inferred_strandness_pattern.txt"
    exit 1
fi

FILE1=$(awk 'BEGIN {FS="\t"} {print $1}' /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/samples.manifest | awk "NR==${SGE_TASK_ID}")
if [[ TRUE == "TRUE" ]]
then
    FILE2=$(awk 'BEGIN {FS="\t"} {print $3}' /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/samples.manifest | awk "NR==${SGE_TASK_ID}")
fi
ID=$(cat /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/samples.manifest | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")
STRANDRULE=$(cat inferred_strandness_pattern.txt)

## Normalizing bigwigs to 40 million 100 bp reads
module load python/2.7.9
module load ucsctools

## Can only use -d when the data is stranded
if [ "${STRANDRULE}" == "none" ]
then
    python ~/.local/bin/bam2wig.py -s /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/hg38.chrom.sizes.gencode -i /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/HISAT2_out/${ID}_accepted_hits.sorted.bam -t 4000000000 -o /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/Coverage/${ID}
elif [ "${STRANDRULE}" == "1++,1--,2+-,2-+" ]
then
    python ~/.local/bin/bam2wig.py -s /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/hg38.chrom.sizes.gencode -i /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/HISAT2_out/${ID}_accepted_hits.sorted.bam -t 4000000000 -o /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/Coverage/${ID} -d "1++,1--,2+-,2-+"
elif [ "${STRANDRULE}" == "1+-,1-+,2++,2--" ]
then
    python ~/.local/bin/bam2wig.py -s /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/hg38.chrom.sizes.gencode -i /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/HISAT2_out/${ID}_accepted_hits.sorted.bam -t 4000000000 -o /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/Coverage/${ID} -d "1+-,1-+,2++,2--"
elif [ "${STRANDRULE}" == "++,--" ]
then
    python ~/.local/bin/bam2wig.py -s /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/hg38.chrom.sizes.gencode -i /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/HISAT2_out/${ID}_accepted_hits.sorted.bam -t 4000000000 -o /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/Coverage/${ID} -d "++,--"
elif [ "${STRANDRULE}" == "+-,-+" ]
then
    python ~/.local/bin/bam2wig.py -s /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/hg38.chrom.sizes.gencode -i /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/HISAT2_out/${ID}_accepted_hits.sorted.bam -t 4000000000 -o /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/Coverage/${ID} -d "+-,-+"
else
    echo "Found unexpected value in inferred_strandness_pattern.txt: ${STRANDRULE}"
fi

## Remove temp files
rm /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/Coverage/${ID}*.wig

echo "**** Job ends ****"
date
