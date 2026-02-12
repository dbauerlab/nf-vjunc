#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process QUANTIFICATION {

    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/r_analysis/${sample}", mode: 'copy', overwrite: true
    container '/nemo/stp/babs/working/bootj/singularity_images/r-analysis.sif'

    input:
        tuple val(sample), path(fastq1), path(fastq2), path(gtf), path(fasta), val(library), path("${sample}.sorted.bam"), path("${sample}.spliced.bam")
        path analysis_script
        path functions_script

    output:
        path "${sample}_*.csv", emit: csv_files
        path "${sample}_*.RDS", emit: rds_files
    
    script:
    """
    Rscript ${analysis_script} ${sample} ${sample}.sorted.bam ${gtf} ${fasta}
    """
    
}