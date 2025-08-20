#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Modules

METADATA {
    take: csv
    main:
        Channel
            .fromPath( csv )
            .splitCsv(header:true)
            .map { row -> [ row.sample, 
                            file(row.fastq1, checkIfExists: true),
                            file(row.fastq2, checkIfExists: true),
                            file(row.gtf, checkIfExists: true),
                            file(row.fasta, checkIfExists: true) ]  }
            .set { data }
    emit:
        data
}

process FLASH {

    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/merged", mode: 'copy', overwrite: true

    container 'staphb/flash:1.2.11'

    input:
        tuple val(sample), path(fastq1), path(fastq2), path(gtf), path(fasta)

    output:
        tuple val(sample), path("${sample}.extendedFrags.fastq.gz"), emit: fastq

    script:
    """
    flash \
      -m 18 \
      -x 0.25 \
      -t ${task.cpus} \
      --allow-outies \
      -z --compress-prog=gzip \
      -o ${sample} \
      -d './'
      ${fastq1} \
      ${fastq2}
    """
}

// Main pipeline

workflow {

    // Run the METADATA workflow
    METADATA(params.input)
    METADATA.out.view { v -> "Channel is ${v}" }

    FLASH(METADATA.out)

}


