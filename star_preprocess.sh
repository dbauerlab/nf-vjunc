#!/usr/bin/env bash

#SBATCH --job-name=star_preprocess
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=ncpu
#SBATCH --time='5:00:00'
#SBATCH --mem=32G
#SBATCH --array=1-19

###############################################################################################################################################
###############################################################################################################################################
## DESCRIPTION:                                                                                                                              ##
## This run of STAR uses a composite host/wuhCor1 genome with standard alignment parameters.                                                 ##
## The input is the trimmed, UMI extracted processed reads collapsed across mates.                                                           ##
## It's function is to allow the creation of a new set of fastq files comprising i) viral and ii) unmapped reads that may be taken forward   ##
## to virus specific alignments steps.                                                                                                       ##
###############################################################################################################################################
###############################################################################################################################################

# Define inputs
THREADS=${SLURM_CPUS_PER_TASK}
PROJDIR=/nemo/stp/babs/working/bootj/projects/bauerd/ziyi.yang/host_virion_rnaseq/richard_pipeline
DESIGN=/nemo/stp/babs/working/bootj/projects/bauerd/ziyi.yang/host_virion_rnaseq/samplesheet.csv
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${DESIGN} | cut -d ',' -f 1)
R1=${FASTQDIR}/$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${DESIGN} | cut -d ',' -f 2)
R2=${FASTQDIR}/$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${DESIGN} | cut -d ',' -f 3)
RESULTS_DIR=${PROJDIR}/star
PREPROCESSED_FASTQ_DIR=${PROJDIR}/fastq_preprocess_outs
VIRAL_FASTQ_DIR=${RESULTS_DIR}/fastq/viral
UNMAPPED_FASTQ_DIR=${RESULTS_DIR}/fastq/unmapped
VIRUS_GENOME=/nemo/stp/babs/working/bootj/genomes/CovBeta/beta.fa
VIRAL_CHR=$(grep '^>' ${VIRUS_GENOME} | sed 's/ .*//' | sed 's/>//' | tr -d '[:space:]' | tr '\n' ' ' | tr '\r' ' ')
STAR_IDX=/nemo/stp/babs/working/bootj/genomes/ChlSab1-1_CovBeta/star

# Make directories 
mkdir -p ${RESULTS_DIR}
mkdir -p ${VIRAL_FASTQ_DIR}
mkdir -p ${UNMAPPED_FASTQ_DIR}

# Load conda module
ml purge
ml Anaconda3/2020.07

# STAR
FQ_FILE=${PREPROCESSED_FASTQ_DIR}/${SAMPLE}.combined.fq.gz

if [ ! -s "${RESULTS_DIR}/${SAMPLE}.Aligned.out.bam" ]
then
  source activate star_2.7.11a
  STAR \
    --runThreadN ${THREADS} \
    --genomeDir ${STAR_IDX} \
    --readFilesIn ${FQ_FILE} \
    --readFilesCommand zcat \
    --twopassMode Basic \
    --outFileNamePrefix ${RESULTS_DIR}/${SAMPLE}. \
    --outReadsUnmapped None \
    --outSAMunmapped Within \
    --outSAMtype BAM Unsorted
  source deactivate
else
  echo "${RESULTS_DIR}/${SAMPLE}.Aligned.out.bam" already exists.
fi

# SAMtools
if [ ! -s "${RESULTS_DIR}/${SAMPLE}.sorted.bam" ]
then
  source activate samtools_1.18
  samtools sort --threads ${THREADS} -o ${RESULTS_DIR}/${SAMPLE}.sorted.bam ${RESULTS_DIR}/${SAMPLE}.Aligned.out.bam
  samtools index ${RESULTS_DIR}/${SAMPLE}.sorted.bam
  samtools idxstats ${RESULTS_DIR}/${SAMPLE}.sorted.bam > ${RESULTS_DIR}/${SAMPLE}.sorted.bam.idxstats
  samtools flagstat ${RESULTS_DIR}/${SAMPLE}.sorted.bam > ${RESULTS_DIR}/${SAMPLE}.sorted.bam.flagstat
  samtools depth -a -m 0 ${RESULTS_DIR}/${SAMPLE}.sorted.bam > ${RESULTS_DIR}/${SAMPLE}.sorted.bam_coverage.txt
  conda deactivate
else
  echo "${RESULTS_DIR}/${SAMPLE}.sorted.bam" already exists.
fi

# Collect reads mapping to viral genome into a new fastq file for further analysis
if [ ! -s "${VIRAL_FASTQ_DIR}/${SAMPLE}.fastq.gz" ]
then
  source activate samtools_1.18  
  samtools view --threads ${THREADS} -b -o ${RESULTS_DIR}/${SAMPLE}.viral.bam ${RESULTS_DIR}/${SAMPLE}.sorted.bam $VIRAL_CHR
  samtools fastq --threads ${THREADS} -0 ${RESULTS_DIR}/${SAMPLE}.viral.fastq ${RESULTS_DIR}/${SAMPLE}.viral.bam
  cat ${RESULTS_DIR}/${SAMPLE}.viral.fastq | gzip > ${VIRAL_FASTQ_DIR}/${SAMPLE}.fastq.gz
  source deactivate
else
  echo "${VIRAL_FASTQ_DIR}/${SAMPLE}.fastq.gz" already exists.
fi

# Collect unmapped reads into a new fastq file for further analysis
if [ ! -s "${UNMAPPED_FASTQ_DIR}/${SAMPLE}.fastq.gz" ]
then
  source activate samtools_1.18  
  samtools view --threads ${THREADS} -b -f 4 -o ${RESULTS_DIR}/${SAMPLE}.unmapped.bam ${RESULTS_DIR}/${SAMPLE}.sorted.bam
  samtools fastq --threads ${THREADS} -0 ${RESULTS_DIR}/${SAMPLE}.unmapped.fastq ${RESULTS_DIR}/${SAMPLE}.unmapped.bam
  cat ${RESULTS_DIR}/${SAMPLE}.unmapped.fastq | gzip > ${UNMAPPED_FASTQ_DIR}/${SAMPLE}.fastq.gz
  source deactivate
else
  echo "${UNMAPPED_FASTQ_DIR}/${SAMPLE}.fastq.gz" already exists.
fi

## Clean up
#if [ -s "${RESULTS_DIR}/${SAMPLE}.sorted.bam" ]
#then 
#  rm ${RESULTS_DIR}/${SAMPLE}.Aligned.out.sam
#  rm ${RESULTS_DIR}/${SAMPLE}.Aligned.out.bam
#fi