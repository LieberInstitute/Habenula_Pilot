#!/bin/bash
#$ -cwd
#$ -l mem_free=32G,h_vmem=34G,h_fsize=100G,bluejay
#$ -N Br1092_Hb-includeintrons
#$ -pe local 4
#$ -o ./logs/Br1092_Hb-includeintrons.txt
#$ -e ./logs/Br1092_Hb-includeintrons.txt
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"
echo "Sample id: Br1092_Hb"
echo "****"


module use /jhpce/shared/jhpce/modulefiles/libd
module load cellranger/6.0.0
cellranger count  --include-introns \
		 --id=Br1092 \
     --transcriptome=/dcl01/ajaffe/data/lab/singleCell/cellranger_reference/refdata-gex-GRCh38-2020-A \
		 --fastqs=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/data_from_linda/2021-05-21/ASpa050421/ \
		 --sample=3c-k \
		 --jobmode=local \
		 --localcores=4 \
		 --localmem=136 \
     