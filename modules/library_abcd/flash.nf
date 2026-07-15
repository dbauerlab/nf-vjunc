#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process FLASH {

    tag "$sample"
    label 'process_medium'

    container 'staphb/flash:1.2.11'

    input:
        tuple val(sample), path(fastq1), path(fastq2)

    output:
        tuple val(sample), path("${sample}.extendedFrags.fastq.gz"), path("${sample}.notCombined_1.fastq.gz"), path("${sample}.notCombined_2.fastq.gz"), emit: mergedfastq

    script:
    """
    flash \
      -m 18 \
      -x 0.25 \
      -t ${task.cpus} \
      --allow-outies \
      -z --compress-prog=gzip \
      --suffix=gz \
      -o ${sample} \
      ${fastq1} \
      ${fastq2}
    """
    
}