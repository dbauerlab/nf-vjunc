#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Concatenates FASTQ files across sequencing lanes for the same sample.
// When only one lane is present, the file is copied as-is.
// gzip multi-stream concatenation (cat file1.gz file2.gz) produces valid gzip output.

process CAT_FASTQ {

    tag "$sample"
    label 'process_low'

    container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"

    input:
        tuple val(sample), path(fastq1s), path(fastq2s), path(gtf), path(fasta), val(library)

    output:
        tuple val(sample), path("${sample}_R1.fastq.gz"), path("${sample}_R2.fastq.gz"), path(gtf), path(fasta), val(library), emit: catfastq

    script:
    """
    cat ${fastq1s} > ${sample}_R1.fastq.gz
    cat ${fastq2s} > ${sample}_R2.fastq.gz
    """

}
