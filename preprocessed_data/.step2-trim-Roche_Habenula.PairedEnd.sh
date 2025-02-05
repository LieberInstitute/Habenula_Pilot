#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=10G,h_vmem=15G,h_fsize=100G
#$ -N step2-trim-Roche_Habenula.PairedEnd
#$ -pe local 3
#$ -o ./logs/trim-Roche_Habenula.$TASK_ID.txt
#$ -e ./logs/trim-Roche_Habenula.$TASK_ID.txt
#$ -t 1-73
#$ -tc 5
#$ -hold_jid pipeline_setup,step1-fastqc-Roche_Habenula.PairedEnd
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

## Locate file and ids
FILE1=$(awk 'BEGIN {FS="\t"} {print $1}' /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/samples.manifest | awk "NR==${SGE_TASK_ID}")
FILEBASE1=$(basename ${FILE1} | sed 's/.fq.gz//; s/.fq//; s/.fastq.gz//; s/.fastq//')
if [ TRUE == "TRUE" ] 
then
    FILE2=$(awk 'BEGIN {FS="\t"} {print $3}' /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/samples.manifest | awk "NR==${SGE_TASK_ID}")
    FILEBASE2=$(basename ${FILE2} | sed 's/.fq.gz//; s/.fq//; s/.fastq.gz//; s/.fastq//')
fi
ID=$(cat /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/samples.manifest | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")

if [ TRUE == "TRUE" ] ; then 
	REPORT1=/dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/FastQC/Untrimmed/${ID}/${FILEBASE1}_fastqc/summary.txt
	REPORT2=/dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/FastQC/Untrimmed/${ID}/${FILEBASE2}_fastqc/summary.txt
	RESULT1=$(grep "Adapter Content" $REPORT1 | cut -c1-4)
	RESULT2=$(grep "Adapter Content" $REPORT2 | cut -c1-4)

	if [[ $RESULT1 == "FAIL" || $RESULT2 == "FAIL" ]] ; then
		## trim, rerun fastQC
		echo "End 1 adapters: $RESULT1"
		echo "End 2 adapters: $RESULT2"
		echo "Trimming will occur."
		
		mkdir -p /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/trimmed_fq
		FP=/dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/trimmed_fq/${ID}_trimmed_forward_paired.fastq
		FU=/dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/trimmed_fq/${ID}_trimmed_forward_unpaired.fastq
		RP=/dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/trimmed_fq/${ID}_trimmed_reverse_paired.fastq
		RU=/dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/trimmed_fq/${ID}_trimmed_reverse_unpaired.fastq
		
		## trim adapters
		java -Xmx5G -jar /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 3 -phred33 		${FILE1} ${FILE2} $FP $FU $RP $RU 		ILLUMINACLIP:/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/Trimmomatic-0.36/adapters/TruSeq2-PE.fa:2:30:10:1 		LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75

		## rerun fastqc
		mkdir -p /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/FastQC/Trimmed/${ID}
		/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/FastQC_v0.11.5/fastqc 		$FP $FU $RP $RU 		--outdir=/dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/FastQC/Trimmed/${ID} --extract
	else
		echo "No trimming required!"
	fi

else
	## reads are single-end
	REPORT1=/dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/FastQC/Untrimmed/${ID}/${FILEBASE1}_fastqc/summary.txt
	RESULT1=$(grep "Adapter Content" $REPORT1 | cut -c1-4)

	if [[ $RESULT1 == "FAIL" ]] ; then
		## trim, rerun fastQC
		echo "Adapters: $RESULT1"
		echo "Trimming will occur."
		
		mkdir -p /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/trimmed_fq
		OUT=/dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/trimmed_fq/${ID}_trimmed.fastq
		
		## trim adapters
		java -Xmx5G -jar /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 3 -phred33 		${FILE1} $OUT 		ILLUMINACLIP:/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/Trimmomatic-0.36/adapters/TruSeq2-SE.fa:2:30:10:1 		LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

		## rerun fastqc
		mkdir -p /dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/FastQC/Trimmed/${ID}
		/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/FastQC_v0.11.5/fastqc $OUT 		--outdir=/dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/FastQC/Trimmed/${ID} --extract
	else
		echo "No trimming required!"
	fi
fi

	
echo "**** Job ends ****"
date
