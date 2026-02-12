#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process STAR_VIRAL {

    tag "$sample"
    label 'process_high'
    publishDir "${params.outdir}/star_viral", mode: 'copy', overwrite: true, pattern: '*.bam'

    container 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_7'

    input:
        tuple val(key), val(sample), path(fastq), path(index)

    output:
        tuple val(sample), path("${sample}_Aligned.out.bam"), emit: viral_bam

    script:
    """
    STAR \
        --runThreadN $task.cpus \
        --genomeDir $index \
        --readFilesIn $fastq \
        --readFilesCommand zcat \
        --outFileNamePrefix ${sample}_ \
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
    """

}