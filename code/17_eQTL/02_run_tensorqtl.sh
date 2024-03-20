#!/bin/bash

#SBATCH -p caracol
#SBATCH -c 1
#SBATCH --gpus=1
#SBATCH --mem=50G
#SBATCH --job-name=02_run_tensorqtl
#SBATCH -o logs/02_run_tensorqtl.log
#SBATCH -e logs/02_run_tensorqtl.log

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load tensorqtl/1.0.8
python 02_run_tensorqtl.py

echo "**** Job ends ****"
date
