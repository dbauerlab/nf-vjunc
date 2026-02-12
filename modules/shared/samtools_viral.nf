#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process SAMTOOLS_VIRAL {
    
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/samtools_viral", mode: 'copy', overwrite: true, pattern: '*.bam'
    publishDir "${params.outdir}/samtools_viral", mode: 'copy', overwrite: true, pattern: '*.fastq.gz'
    publishDir "${params.outdir}/samtools_viral", mode: 'copy', overwrite: true, pattern: '*.idxstats'
    publishDir "${params.outdir}/samtools_viral", mode: 'copy', overwrite: true, pattern: '*.flagstat'
    publishDir "${params.outdir}/samtools_viral", mode: 'copy', overwrite: true, pattern: '*.txt'

    container 'quay.io/biocontainers/samtools:1.22.1--h96c455f_0'

    input:
        tuple val(sample), path(bam)

    output:
        tuple val(sample), path("${sample}.idxstats"), path("${sample}.flagstat"), path("${sample}.coverage.txt"), emit: stats
        tuple val(sample), path("${sample}.sorted.bam"), path("${sample}.spliced.bam"), emit: bams

    script:
    """
    # Sort, index, and generate stats on the pre-mapping BAM file
    samtools sort --threads $task.cpus -o ${sample}.sorted.bam $bam
    samtools index ${sample}.sorted.bam
    samtools idxstats ${sample}.sorted.bam > ${sample}.idxstats
    samtools flagstat ${sample}.sorted.bam > ${sample}.flagstat
    samtools depth -a -m 0 ${sample}.sorted.bam > ${sample}.coverage.txt

    # Collect reads that are spliced
    samtools view -h ${sample}.sorted.bam | awk -v OFS="\t" '\$0 ~ /^@/{print \$0;next;} \$6 ~ /N/' | samtools view -b -o ${sample}.spliced.bam
    samtools index ${sample}.spliced.bam
    """
    
}