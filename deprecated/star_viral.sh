#!/usr/bin/env bash

#SBATCH --job-name=STAR
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=ncpu
#SBATCH --time='72:00:00'
#SBATCH --mem=256G
#SBATCH --array=1-19

# If jobs fail due to memory constraints, then change the following and rerun:
# --cpus-per-task=32
# --mem=250G

# Define inputs
THREADS=${SLURM_CPUS_PER_TASK}
PROJDIR=/nemo/stp/babs/working/bootj/projects/bauerd/ziyi.yang/host_virion_rnaseq/richard_pipeline
DESIGN=/nemo/stp/babs/working/bootj/projects/bauerd/ziyi.yang/host_virion_rnaseq/samplesheet.csv
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${DESIGN} | cut -d ',' -f 1)
RESULTS_DIR=${PROJDIR}/star/viral_results
VIRAL_FASTQ_DIR=${PROJDIR}/star/fastq/viral
STAR_IDX=/nemo/stp/babs/working/bootj/genomes/CovBeta/star
STAR_GENOME=/nemo/stp/babs/working/bootj/genomes/CovBeta/beta.fa

# Directories
mkdir -p ${RESULTS_DIR}

# Load conda module
ml purge
ml Anaconda3/2020.07

# STAR
FQ_FILE="${VIRAL_FASTQ_DIR}/${SAMPLE}.fastq.gz"

if [ ! -s "${RESULTS_DIR}${SAMPLE}.Aligned.out.bam" ]
then
  source activate star_2.7.11a
  STAR \
    --runThreadN ${THREADS} \
    --genomeDir ${STAR_IDX} \
    --readFilesIn ${FQ_FILE} \
    --readFilesCommand zcat \
    --outFileNamePrefix ${RESULTS_DIR}/${SAMPLE}. \
    --outReadsUnmapped Fastx \
    --outStd Log \
    --outSAMtype BAM Unsorted \
    --outSAMattributes Standard \
    --twopassMode Basic \
    --seedPerWindowNmax 30 \
    --alignIntronMin 1 \
    --outSJfilterOverhangMin 20 20 20 20 \
    --outSJfilterCountUniqueMin 1 1 1 1 \
    --outSJfilterCountTotalMin 1 1 1 1 \
    --outSJfilterDistToOtherSJmin 0 0 0 0 \
    --scoreGapNoncan 0 \
    --scoreGapGCAG 0 \
    --scoreGapATAC 0 \
    --alignSJoverhangMin 20 \
    --outFilterMatchNmin 40 \
    --outSJfilterReads All \
    --outSAMmultNmax 1 \
    --outFilterMismatchNoverLmax 0.1 \
    --alignEndsType Local \
    --outFilterType BySJout \
    --limitOutSJcollapsed 10000000 \
    --limitIObufferSize 1500000000 1500000000 \
    --alignSJstitchMismatchNmax 0 0 0 0 \
    --alignSJDBoverhangMin 20 \
    --alignSoftClipAtReferenceEnds Yes \
    --scoreGenomicLengthLog2scale 0 \
    --outFilterMultimapNmax 1

  conda deactivate
else
  echo "${RESULTS_DIR}/${SAMPLE}.Aligned.out.bam" already exists.
fi


# SAMtools
if [ ! -s "${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.bam" ]
then
  source activate samtools_1.18
  samtools sort --threads $THREADS -o ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.bam ${RESULTS_DIR}/${SAMPLE}.Aligned.out.bam
  samtools index ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.bam
  samtools idxstats ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.bam > ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.bam.idxstats
  samtools flagstat ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.bam > ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.bam.flagstat
  samtools depth -a -m 0 ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.bam > ${RESULTS_DIR}/${SAMPLE}.sorted.bam_coverage.txt
  conda deactivate
else
  echo "${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.bam" already exists.
fi

# SAMtools, subset spliced
if [ ! -s "${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.spliced.bam" ]
then
  source activate samtools_1.18
  samtools view -h ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.bam | awk -v OFS="\t" '$0 ~ /^@/{print $0;next;} $6 ~ /N/' | samtools view -b -o ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.spliced.bam
  samtools index ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.spliced.bam
  conda deactivate
else
  echo "${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.spliced.bam" already exists.
fi

source activate samtools_1.18
CPM_FACTOR_ALLREADS=$(bc <<< "scale=6;1000000/$(samtools view -f 0 -c ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.bam)")
CPM_FACTOR_SPLICEDREADS=$(bc <<< "scale=6;1000000/$(samtools view -f 0 -c ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.spliced.bam)")
conda deactivate


## BEDtools
if [ ! -s "${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.unstranded.bigwig" ]
then
  source activate bedtools_2.31.1
  bedtools genomecov \
    -scale $CPM_FACTOR_ALLREADS \
    -bga \
    -ibam ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.bam > ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.unstranded.bg
  ${PROJDIR}/scripts/bedSort ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.unstranded.bg ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.unstranded.sorted.bg
  ${PROJDIR}/scripts/bedGraphToBigWig ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.unstranded.sorted.bg ${STAR_GENOME}.fai ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.unstranded.bigwig
  bedtools genomecov \
    -scale $CPM_FACTOR_ALLREADS \
    -bga \
    -strand "+" \
    -ibam ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.bam > ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.plus.bg
  ${PROJDIR}/scripts/bedSort ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.plus.bg ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.plus.sorted.bg
  ${PROJDIR}/scripts/bedGraphToBigWig ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.plus.sorted.bg ${STAR_GENOME}.fai ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.plus.bigwig
  bedtools genomecov \
    -scale $CPM_FACTOR_ALLREADS \
    -bga \
    -strand "-" \
    -ibam ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.bam > ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.minus.bg
  ${PROJDIR}/scripts/bedSort ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.minus.bg ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.minus.sorted.bg
  ${PROJDIR}/scripts/bedGraphToBigWig ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.minus.sorted.bg ${STAR_GENOME}.fai ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.minus.bigwig

  rm ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.unstranded.bg
  rm ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.plus.bg
  rm ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.plus.sorted.bg
  rm ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.minus.bg
  rm ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.unstranded.sorted.bg
  rm ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.minus.sorted.bg

  conda deactivate
else
  echo "${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.unstranded.bigwig" already exists.
fi


## BEDtools
if [ ! -s "${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.spliced.unstranded.bigwig" ]
then
 source activate bedtools_2.31.1
  bedtools genomecov \
    -scale $CPM_FACTOR_ALLREADS \
    -bga \
    -ibam ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.spliced.bam > ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.spliced.unstranded.bg
  ${PROJDIR}/scripts/bedSort ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.spliced.unstranded.bg ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.spliced.unstranded.sorted.bg
  ${PROJDIR}/scripts/bedGraphToBigWig ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.spliced.unstranded.sorted.bg ${STAR_GENOME}.fai ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.spliced.unstranded.bigwig
  bedtools genomecov \
    -scale $CPM_FACTOR_ALLREADS \
    -bga \
    -strand "+" \
    -ibam ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.spliced.bam > ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.spliced.plus.bg
  ${PROJDIR}/scripts/bedSort ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.spliced.plus.bg ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.spliced.plus.sorted.bg
  ${PROJDIR}/scripts/bedGraphToBigWig ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.spliced.plus.sorted.bg ${STAR_GENOME}.fai ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.spliced.plus.bigwig
  bedtools genomecov \
    -scale $CPM_FACTOR_ALLREADS \
    -bga \
    -strand "-" \
    -ibam ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.spliced.bam > ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.spliced.minus.bg
  ${PROJDIR}/scripts/bedSort ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.spliced.minus.bg ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.spliced.minus.sorted.bg
  ${PROJDIR}/scripts/bedGraphToBigWig ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.spliced.minus.sorted.bg ${STAR_GENOME}.fai ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.spliced.minus.bigwig

  rm ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.spliced.unstranded.bg
  rm ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.spliced.plus.bg
  rm ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.spliced.minus.bg
  rm ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.spliced.unstranded.sorted.bg
  rm ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.spliced.plus.sorted.bg
  rm ${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.spliced.minus.sorted.bg

  conda deactivate
else
  echo "${RESULTS_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.spliced.unstranded.bigwig" already exists.
fi