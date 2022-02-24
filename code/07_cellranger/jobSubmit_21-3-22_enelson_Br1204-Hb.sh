#!/bin/bash
#$ -cwd
#$ -l mem_free=32G,h_vmem=34G,h_fsize=100G,bluejay
#$ -N test_Br2_Hb-includeintrons
#$ -pe local 4
#$ -o ./logs/Br1204_Hb-includeintrons.txt
#$ -e ./logs/Br1204_Hb-includeintrons.txt
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"
echo "Sample id: Br1204_Hb"
echo "****"


module use /jhpce/shared/jhpce/modulefiles/libd
module load cellranger/6.0.0
cellranger count  --include-introns \
		 --id=Br1204 \
     --transcriptome=/dcl01/ajaffe/data/lab/singleCell/cellranger_reference/refdata-gex-GRCh38-2020-A \
		 --fastqs=/dcl01/ajaffe/data/lab/singleCell/10x_pilot/FASTQ/Feb2021/Br1204_Hb/ \
		 --sample=Br1204-Hb \
		 --jobmode=local \
		 --localcores=4 \
		 --localmem=136 \