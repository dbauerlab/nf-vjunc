#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process SAMTOOLS_HOST {
    
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/samtools_host", mode: 'copy', overwrite: true, pattern: '*.bam'
    publishDir "${params.outdir}/samtools_host", mode: 'copy', overwrite: true, pattern: '*.fastq.gz'
    publishDir "${params.outdir}/samtools_host", mode: 'copy', overwrite: true, pattern: '*.idxstats'
    publishDir "${params.outdir}/samtools_host", mode: 'copy', overwrite: true, pattern: '*.flagstat'
    publishDir "${params.outdir}/samtools_host", mode: 'copy', overwrite: true, pattern: '*.coverage.txt'

    container 'quay.io/biocontainers/samtools:1.22.1--h96c455f_0'

    input:
        tuple val(sample), path(fasta), path(gtf), path(bam)

    output:
        tuple val(sample), path("${sample}.idxstats"), path("${sample}.flagstat"), path("${sample}.coverage.txt"), emit: stats
        tuple val(sample), path(fasta), path(gtf), path("${sample}.viral.bam"), path("${sample}.viral.fastq.gz"), emit: viral
        tuple val(sample), path("${sample}.unmapped.bam"), path("${sample}.unmapped.fastq.gz"), emit: unmapped

    script:
    """
    VIRAL_CHR=\$(grep '^>' ${fasta} | sed 's/ .*//' | sed 's/>//' | tr -d '[:space:]' | tr '\n' ' ' | tr '\r' ' ')

    # Sort, index, and generate stats on the pre-mapping BAM file
    samtools sort --threads $task.cpus -o ${sample}.sorted.bam $bam
    samtools index ${sample}.sorted.bam
    samtools idxstats ${sample}.sorted.bam > ${sample}.idxstats
    samtools flagstat ${sample}.sorted.bam > ${sample}.flagstat
    samtools depth -a -m 0 ${sample}.sorted.bam > ${sample}.coverage.txt

    # Collect reads mapping to viral genome into a new fastq file for further analysis
    samtools view --threads $task.cpus -b -o ${sample}.viral.bam ${sample}.sorted.bam \${VIRAL_CHR}
    samtools fastq --threads $task.cpus -0 ${sample}.viral.fastq.gz ${sample}.viral.bam

    # Collect unmapped reads into a new fastq file for further analysis
    samtools view --threads $task.cpus -b -f 4 -o ${sample}.unmapped.bam ${sample}.sorted.bam
    samtools fastq --threads $task.cpus -0 ${sample}.unmapped.fastq.gz ${sample}.unmapped.bam
    """
    
}