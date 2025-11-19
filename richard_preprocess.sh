#!/usr/bin/env bash
#
#SBATCH --job-name=preprocess_reads
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=ncpu
#SBATCH --time='72:00:00'
#SBATCH --mem=32G
#SBATCH --array=1-60

## This script takes the gathered fastq files (i.e. samples concatenated across lanes) and performs library-specific pre-processing:
# Adapter trimming
# UMI extraction
# Hard-clipping
# Collapse across mates

THREADS=${SLURM_CPUS_PER_TASK}

DESIGN="/nemo/stp/babs/working/mitterr/projects/bauerd/nuno.santos/PM23265-ns578/data/experiment_table.csv"

IFS=',' read LIMSID VIRUS HPI REP LIB DATE LABEL HOST NOTES FQ1 FQ2 <<< $(awk --assign LN=$((${SLURM_ARRAY_TASK_ID}+1)) 'NR==LN{print}' $DESIGN)


## Params
ADAPTER="AGATCGGAAGAGC"
UMI_R1="NNNNNNNNNN"
UMI_R2="NNNNNNNNNN"
LIB_B_PCRPRIMER_HARDCLIP=27
LIB_C_PCRPRIMER_HARDCLIP=19
LIB_D_PCRPRIMER_HARDCLIP=19
MIN_OVERLAP=18


## Directories
PROJDIR="/nemo/stp/babs/working/mitterr/projects/bauerd/nuno.santos/PM23265-ns578/"
SCRATCHDIR="/flask/scratch/babs/mitterr/projects/bauerd/nuno.santos/PM23265-ns578/"
RESULTSDIR="${PROJDIR}results/"
TMPDIR="${SCRATCHDIR}tmp/${LIMSID}/"
#FASTQDIR="${RESULTSDIR}fastq/preprocessed/"
FASTQDIR="${SCRATCHDIR}fastq/preprocessed/"
mkdir -p $PROJDIR
mkdir -p $TMPDIR
mkdir -p $RESULTSDIR
mkdir -p $FASTQDIR


## Load conda module
ml purge
ml Anaconda3/2020.07


## 3’ adapter trimming using AGATCGGAAGAGC for both R1 and R2
if [ ! -s "${TMPDIR}${LIMSID}_val_1.fq.gz" ]
then
  source activate trim-galore_0.6.10
  trim_galore \
    --paired \
    -a ${ADAPTER} \
    -a2 ${ADAPTER} \
    --cores ${THREADS} \
    --output_dir ${TMPDIR} \
    --basename ${LIMSID} \
    --length 30 \
    --fastqc \
    ${FQ1} \
    ${FQ2}
  conda deactivate
else
  echo ${TMPDIR}${LIMSID}_val_1.fq.gz already exists.
fi


## Extract UMI as the first 10bp of R1 and of R2.  Extracted sequence will be moved to the FASTQ headers and hard-clipped from the actual sequence.
if [ ! -s "${TMPDIR}${LIMSID}.umi_extracted.R1.fastq.gz" ]
then
  source activate umi_tools_1.1.4
  umi_tools extract \
   -I ${TMPDIR}${LIMSID}_val_1.fq.gz \
   --read2-in=${TMPDIR}${LIMSID}_val_2.fq.gz \
   --bc-pattern=${UMI_R1} \
   --bc-pattern2=${UMI_R2} \
   --stdout=${TMPDIR}${LIMSID}.umi_extracted.R1.fastq.gz \
   --read2-out=${TMPDIR}${LIMSID}.umi_extracted.R2.fastq.gz \
   --log=${TMPDIR}${LIMSID}.umi_log
  conda deactivate
else
  echo ${TMPDIR}${LIMSID}.umi_extracted.R1.fastq.gz already exits.
fi


## Hard-clip the PCR primer from the sequences downstream of the removed R1 UMI.  The amount of the sequence to clip is library specific.  For library B the upper limit of the expected size is 27bp.
if [ ! -s "${TMPDIR}${LIMSID}.to_collapse.R1.fq.gz" ]
then
  if test "${LIB}" != "A"
  then 
    HARDCLIP=0
    if test "${LIB}" = "B"
    then
      HARDCLIP=$LIB_B_PCRPRIMER_HARDCLIP
    fi
    if test "${LIB}" = "C"
    then
      HARDCLIP=$LIB_C_PCRPRIMER_HARDCLIP
    fi
    if test "${LIB}" = "D"
    then
      HARDCLIP=$LIB_D_PCRPRIMER_HARDCLIP
    fi
    HARDCLIP=$(( $HARDCLIP + 1 ))
    source activate fastx_toolkit_0.0.14
    zcat ${TMPDIR}${LIMSID}.umi_extracted.R1.fastq.gz | fastx_trimmer \
      -z \
      -f $HARDCLIP \
      -o ${TMPDIR}${LIMSID}.umi_extracted.clipped.R1.gz
    conda deactivate

    ln -s ${TMPDIR}${LIMSID}.umi_extracted.clipped.R1.gz ${TMPDIR}${LIMSID}.to_collapse.R1.fq.gz
    ln -s ${TMPDIR}${LIMSID}.umi_extracted.R2.fastq.gz ${TMPDIR}${LIMSID}.to_collapse.R2.fq.gz
  else
    ln -s ${TMPDIR}${LIMSID}.umi_extracted.R1.fastq.gz ${TMPDIR}${LIMSID}.to_collapse.R1.fq.gz
    ln -s ${TMPDIR}${LIMSID}.umi_extracted.R2.fastq.gz ${TMPDIR}${LIMSID}.to_collapse.R2.fq.gz
  fi
else
  echo ${TMPDIR}${LIMSID}.to_collapse.R1.fq.gz already exists.
fi


## Collapse overlapping R1 and R2 reads based on a defined minimum overlap.
if [ ! -s "${TMPDIR}${LIMSID}.extendedFrags.fastq.gz" ]
then
  source activate flash_1.2.11
    flash \
      -m 18 \
      -x 0.25 \
      -t $THREADS \
      --allow-outies \
      -z --compress-prog=gzip \
      -o ${LIMSID} \
      -d ${TMPDIR} \
      ${TMPDIR}${LIMSID}.to_collapse.R1.fq.gz \
      ${TMPDIR}${LIMSID}.to_collapse.R2.fq.gz
  conda deactivate
else
  echo ${TMPDIR}${LIMSID}.extendedFrags.fastq.gz already exists.
fi


## Combine R1R2 collapsed reads with R1 singletones (LIB B or C or D), or R1R2 + R1 + R2 (LIB A)
## R2 singletons are first reverse-complement to put them in the same orientation as R1 singletones.
FQ_COMBINED="${TMPDIR}${LIMSID}.combined.fq.gz"
if [ ! -s ${FQ_COMBINED} ]
then
  if test "${LIB}" == "A"
  then
    source activate fastx_toolkit_0.0.14
    zcat ${TMPDIR}${LIMSID}.notCombined_2.fastq.gzip | fastx_reverse_complement \
      -z \
      -o ${TMPDIR}${LIMSID}.notCombined_2.fastq.revcomp.gzip
    conda deactivate

    cat ${TMPDIR}${LIMSID}.extendedFrags.fastq.gzip ${TMPDIR}${LIMSID}.notCombined_1.fastq.gzip ${TMPDIR}${LIMSID}.notCombined_2.fastq.revcomp.gzip > $FQ_COMBINED
  else
    cat ${TMPDIR}${LIMSID}.extendedFrags.fastq.gzip ${TMPDIR}${LIMSID}.notCombined_1.fastq.gzip > $FQ_COMBINED
  fi 
else
  echo $FQ_COMBINED already exists.
fi


## Reverse complement
FQ_COMBINED_REVCOMP="${TMPDIR}${LIMSID}.combined.revcomp.fq.gz"
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


## Rearrange the fastq headers to keep the UMI at the end separated by an underscore.  The rest is placed before the UMI, separated from the rest of the read name by a backslash.
## This is to prevent trimming of the headers resulting in dupicate names from R1 and R2 pairs during subsequent alignment.
FQ_FILE=${TMPDIR}${LIMSID}.combined.revcomp.adjusted_header.fq.gz
if [ ! -s $FQ_FILE ]
then
  zcat $FQ_COMBINED_REVCOMP | sed --regexp-extended 's/(^@\S+)_(\S+) ((.):\S+$)/\1\\\3_\2/' | gzip > $FQ_FILE
fi


## Remove duplicates entries based on read headers.  This is an artifact of 


### Copy processed fastq files
FQ_FILE=${FASTQDIR}${LIMSID}.fastq.gz
if [ ! -s $FQ_FILE ]
then
  cp ${TMPDIR}${LIMSID}.combined.revcomp.adjusted_header.fq.gz $FQ_FILE
fi
