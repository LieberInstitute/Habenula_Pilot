#!/bin/bash

## Usage:
# sh 01_get_droplet_scores_loop.sh

## Create the logs directory
mkdir -p logs

for sample in Br1092 Br1204 Br1469 Br1735 Br5555 Br5558 Br5639; do

    ## Internal script name
    SHORT="01_get_droplet_scores_loop_${sample}"
    NAME="get_droplet_scores_loop_${sample}"

    # Construct shell file
    echo "Creating script 01_get_droplet_scores_loop_${sample}"
    cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=50G,h_vmem=50G,h_fsize=100G
#$ -N ${NAME}
#$ -o logs/${SHORT}.txt
#$ -e logs/${SHORT}.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "Task id: \${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 01_get_droplet_scores.R ${sample}

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/


EOF

    call="qsub .${SHORT}.sh"
    echo $call
    $call
done
