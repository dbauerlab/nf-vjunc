#!/usr/bin/env bash

#SBATCH --job-name=preprocess_reads
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=ncpu
#SBATCH --time='72:00:00'
#SBATCH --mem=32G
#SBATCH --array=1-19

# This script builds on Richard Mitter's original script for Ziyi's Host RNA virion packaging project
# Takes the gathered fastq files and performs library-specific pre-processing
# Trimming already performed by George Young - use these files as input
# Then just: 
# Collapse across mates

# Define inputs
THREADS=${SLURM_CPUS_PER_TASK}
PROJDIR=/nemo/stp/babs/working/bootj/projects/bauerd/ziyi.yang/host_virion_rnaseq
DESIGN=${PROJDIR}/samplesheet.csv
FASTQDIR=/nemo/lab/bauerd/data/STPs/babs/inputs/ziyi.yang/george.young/220831_Bauer_CoV-virion-RNA/trimmed
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${DESIGN} | cut -d ',' -f 1)
R1=${FASTQDIR}/$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${DESIGN} | cut -d ',' -f 2)
R2=${FASTQDIR}/$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${DESIGN} | cut -d ',' -f 3)
RESULTSDIR=${PROJDIR}/richard_pipeline/fastq_preprocess_outs

# Create directories
mkdir -p ${RESULTSDIR}

# Load conda module
ml purge
ml Anaconda3/2020.07

# Collapse overlapping R1 and R2 reads based on a defined minimum overlap.
if [ ! -s "${RESULTSDIR}/${SAMPLE}.extendedFrags.fastq.gz" ]
then
  source activate flash_1.2.11
    flash \
      -m 18 \
      -x 0.25 \
      -t ${THREADS} \
      --allow-outies \
      -z --compress-prog=gzip \
      -o ${SAMPLE} \
      -d ${RESULTSDIR} \
      ${R1} \
      ${R2}
  conda deactivate
else
  echo ${RESULTSDIR}/${SAMPLE}.extendedFrags.fastq.gz already exists.
fi

# Combine R1+R2 collapsed reads with R1 singletons (collapsed read is on R1 strand so no need to edit R1 singletons)
# Also combine R1+R2 collapsed reads with R2 singletons (R2 singletons need to be reverse complemented)
FQ_COMBINED=${RESULTSDIR}/${SAMPLE}.combined.fq.gz
if [ ! -s ${FQ_COMBINED} ]
then
    source activate fastx_toolkit_0.0.14
    zcat ${RESULTSDIR}/${SAMPLE}.notCombined_2.fastq.gzip | fastx_reverse_complement \
      -z \
      -o ${RESULTSDIR}/${SAMPLE}.notCombined_2.fastq.revcomp.gzip
    conda deactivate

    cat ${RESULTSDIR}/${SAMPLE}.extendedFrags.fastq.gzip ${RESULTSDIR}/${SAMPLE}.notCombined_1.fastq.gzip ${RESULTSDIR}/${SAMPLE}.notCombined_2.fastq.revcomp.gzip > $FQ_COMBINED
else
  echo ${FQ_COMBINED} already exists.
fi

## Reverse complement
FQ_COMBINED_REVCOMP="${RESULTSDIR}/${SAMPLE}.combined.revcomp.fq.gz"
if [ ! -s "${FQ_COMBINED_REVCOMP}" ]
then
  source activate fastx_toolkit_0.0.14
  zcat $FQ_COMBINED | fastx_reverse_complement \
    -z \
    -o $FQ_COMBINED_REVCOMP
  conda deactivate
else
  echo ${FQ_COMBINED_REVCOMP} already exists.
fi