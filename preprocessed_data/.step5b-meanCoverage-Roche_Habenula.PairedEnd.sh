#!/bin/bash
#$ -cwd
#$ -l leek,mem_free=100G,h_vmem=120G,h_fsize=100G
#$ -N step5b-meanCoverage-Roche_Habenula.PairedEnd
#$ -o ./logs/meanCoverage-Roche_Habenula.txt
#$ -e ./logs/meanCoverage-Roche_Habenula.txt
#$ -hold_jid pipeline_setup,step5-coverage-Roche_Habenula.PairedEnd
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

if [ ! -f "inferred_strandness_pattern.txt" ]
then
    echo "Missing the file inferred_strandness_pattern.txt"
    exit 1
fi

## Force R 3.3.x in JHPCE (to avoid some issues with conda_R)
module unload conda_R
module load R/3.3.x

## Load required software
module load wiggletools/default
module load ucsctools

## Read the strandness information
STRANDRULE=$(cat inferred_strandness_pattern.txt)

if [ ${STRANDRULE} == "none" ]
then
    echo "*****************************"
    date
    echo "Processing unstranded data"
    echo "*****************************"
    ## Locate normalized BigWig files and concatenate them in a space separated list
    BIGWIGS=$(cat /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/samples.manifest | awk '{print "/dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/Coverage/"$NF".bw"}' | paste -sd " ")
    
    ## Create mean of normalized bigwigs
    wiggletools write /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/Coverage/mean.wig mean ${BIGWIGS}
    wigToBigWig /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/Coverage/mean.wig /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/hg38.chrom.sizes.gencode /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/Coverage/mean.bw
else
    for strand in Forward Reverse
    do
        echo "*****************************"
        date
        echo "Processing strand ${strand}"
        echo "*****************************"
        ## Locate normalized BigWig files and concatenate them in a space separated list
        BIGWIGS=$(cat /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/samples.manifest | awk -v strand="${strand}" '{print "/dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/Coverage/"$NF"."strand".bw"}' | paste -sd " ")
    
        ## Create mean of normalized bigwigs
        wiggletools write /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/Coverage/mean.${strand}.wig mean ${BIGWIGS}
        wigToBigWig /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/Coverage/mean.${strand}.wig /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/hg38.chrom.sizes.gencode /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/Coverage/mean.${strand}.bw
    done
fi

## Remove temp files
rm /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/Coverage/mean*.wig

echo "**** Job ends ****"
date
