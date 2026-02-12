#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process STAR_HOST {

    tag "$sample"
    label 'process_high'
    publishDir "${params.outdir}/star_host", mode: 'copy', overwrite: true, pattern: '*.bam'

    container 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_7'

    input:
        tuple val(key), val(sample), path(combined), path(reverse), path(fasta), path(gtf), path(joint_index)

    output:
        tuple val(sample), path(fasta), path(gtf), path("${sample}_Aligned.out.bam"), emit: host_bam

    script:
    """
    STAR \
        --runThreadN $task.cpus \
        --genomeDir $joint_index \
        --readFilesIn $reverse \
        --readFilesCommand zcat \
        --twopassMode Basic \
        --outReadsUnmapped None \
        --outSAMunmapped Within \
        --outSAMtype BAM Unsorted \
        --outFileNamePrefix ${sample}_
    """

}