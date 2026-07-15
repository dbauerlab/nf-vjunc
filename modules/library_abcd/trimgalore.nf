#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process TRIMGALORE {

    tag "$sample"
    label 'process_medium'

    container 'quay.io/biocontainers/trim-galore:0.6.9--hdfd78af_0'

    input:
        tuple val(sample), path(fastq1), path(fastq2), path(gtf), path(fasta), val(library)

    output:
        tuple val(sample), path("${sample}_val_1.fq.gz"), path("${sample}_val_2.fq.gz"), emit: trimfastq

    script:
    """
    trim_galore \
        --paired \
        -a "AGATCGGAAGAGC" \
        -a2 "AGATCGGAAGAGC" \
        --cores ${task.cpus} \
        --basename ${sample} \
        --length 30 \
        ${fastq1} \
        ${fastq2}
    """
    
}