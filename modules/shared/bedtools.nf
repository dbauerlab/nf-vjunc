#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process BEDTOOLS {
    
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/bedtools", mode: 'copy', overwrite: true, pattern: '*.sorted.bg'

    container 'quay.io/biocontainers/bedtools:2.31.1--h13024bc_3'

    input:
        tuple val(sample), path(sorted), path(spliced)

    output:
        tuple val(sample), path("${sample}.unstranded.cov.sorted.bg"), path("${sample}.plus.cov.sorted.bg"), path("${sample}.minus.cov.sorted.bg"), path("${sample}.spliced.unstranded.cov.sorted.bg"), path("${sample}.spliced.plus.cov.sorted.bg"), path("${sample}.spliced.minus.cov.sorted.bg"), emit: bedgraph
        
    script:
    """
    # Do coverage for reads from unfiltered BAM
    bedtools genomecov \
        -bga \
        -ibam ${sorted} > ${sample}.unstranded.cov.bg
    bedtools sort -i ${sample}.unstranded.cov.bg > ${sample}.unstranded.cov.sorted.bg
  
    bedtools genomecov \
        -bga \
        -strand "+" \
        -ibam ${sorted} > ${sample}.plus.cov.bg
    bedtools sort -i ${sample}.plus.cov.bg > ${sample}.plus.cov.sorted.bg

    bedtools genomecov \
        -bga \
        -strand "-" \
        -ibam ${sorted} > ${sample}.minus.cov.bg
    bedtools sort -i ${sample}.minus.cov.bg > ${sample}.minus.cov.sorted.bg
  
    # Do coverage for reads from spliced BAM
    bedtools genomecov \
        -bga \
        -ibam ${spliced} > ${sample}.spliced.unstranded.cov.bg
    bedtools sort -i ${sample}.spliced.unstranded.cov.bg > ${sample}.spliced.unstranded.cov.sorted.bg
  
    bedtools genomecov \
        -bga \
        -strand "+" \
        -ibam ${spliced} > ${sample}.spliced.plus.cov.bg
    bedtools sort -i ${sample}.spliced.plus.cov.bg > ${sample}.spliced.plus.cov.sorted.bg

    bedtools genomecov \
        -bga \
        -strand "-" \
        -ibam ${spliced} > ${sample}.spliced.minus.cov.bg
    bedtools sort -i ${sample}.spliced.minus.cov.bg > ${sample}.spliced.minus.cov.sorted.bg
    """
    
}