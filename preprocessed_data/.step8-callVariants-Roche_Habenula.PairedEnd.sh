#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=2G,h_vmem=3G,h_fsize=100G
#$ -N step8-callVariants-Roche_Habenula.PairedEnd
#$ -o ./logs/callVariants-Roche_Habenula.$TASK_ID.txt
#$ -e ./logs/callVariants-Roche_Habenula.$TASK_ID.txt
#$ -t 7
#$ -tc 100
#$ -hold_jid pipeline_setup,step3-hisat2-Roche_Habenula.PairedEnd
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

mkdir -p /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/Genotypes/

ID=$(cat /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/samples.manifest | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")
BAM=/dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/HISAT2_out/${ID}_accepted_hits.sorted.bam

SNPTMP=/dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/Genotypes/${ID}_tmp.vcf
SNPOUTGZ=/dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/Genotypes/${ID}.vcf.gz

module load bcftools
/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/samtools-1.2/samtools mpileup -l /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Genotyping/common_missense_SNVs_hg38.bed -AB -q0 -Q13 -d1000000 -uf /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/GRCh38.primary_assembly.genome.fa ${BAM} -o ${SNPTMP}
bcftools call -mv -Oz ${SNPTMP} > ${SNPOUTGZ}

module load htslib
tabix -p vcf ${SNPOUTGZ}

rm ${SNPTMP}

echo "**** Job ends ****"
date
