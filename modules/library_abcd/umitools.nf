#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process UMITOOLS {

    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/umitools", mode: 'copy', overwrite: true

    container 'quay.io/biocontainers/umi_tools:1.1.5--py39hbcbf7aa_4'

    input:
        tuple val(sample), path(fastq1), path(fastq2)

    output:
        tuple val(sample), path("${sample}.umi_extracted.R1.fastq.gz"), path("${sample}.umi_extracted.R2.fastq.gz"), emit: umifastq
        tuple val(sample), path("${sample}.umi_log"), emit: umilog

    script:
    """
    umi_tools extract \
        -I ${sample}_val_1.fq.gz \
        --read2-in=${sample}_val_2.fq.gz \
        --bc-pattern="NNNNNNNNNN" \
        --bc-pattern2="NNNNNNNNNN" \
        --stdout=${sample}.umi_extracted.R1.fastq.gz \
        --read2-out=${sample}.umi_extracted.R2.fastq.gz \
        --log=${sample}.umi_log
    """
    
}