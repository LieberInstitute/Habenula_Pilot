#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=20G
#SBATCH --job-name=03_adj_p_vals
#SBATCH -o /dev/null
#SBATCH -e /dev/null

run_mode=interaction
interaction_cov=tot_Hb

if [[ "$run_mode" == "interaction" ]]; then
    log_path="logs/03_adj_p_val_${run_mode}_${interaction_cov}.log"
else
    log_path="logs/03_adj_p_val_${run_mode}.log"
fi

{
set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load conda_R/4.3.x
Rscript 03_adj_p_vals.R -m $run_mode -c $interaction_cov

echo "**** Job ends ****"
date
} > $log_path 2>&1
