#!/bin/bash -l
#SBATCH --output=logs/05_Hashikawa_export_markers.txt
#SBATCH --error=logs/05_Hashikawa_export_markers.txt
#SBATCH --partition=shared
#SBATCH --job-name=0Hashikawa_export_markers
#SBATCH --mem=5GB

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## load modules
module load conda_R/4.3
# list modules
module list

Rscript 05_Hashikawa_export_markers.R

echo "**** Job ends ****"
date