#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=8G,h_vmem=10G
#$ -N step8b-mergeVariantCalls-Roche_Habenula.PairedEnd
#$ -o ./logs/mergeVariantCalls-Roche_Habenula.txt
#$ -e ./logs/mergeVariantCalls-Roche_Habenula.txt
#$ -hold_jid pipeline_setup,step8-callVariants-Roche_Habenula.PairedEnd
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

VCFS=$(cat /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/samples.manifest | awk '{print "/dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/Genotypes/"$NF".vcf.gz"}' | paste -sd " ")

module load vcftools
module load htslib
vcf-merge ${VCFS} | bgzip -c > /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/Genotypes/mergedVariants.vcf.gz
tabix -p vcf /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/Genotypes/mergedVariants.vcf.gz

echo "**** Job ends ****"
date
