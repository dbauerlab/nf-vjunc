#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Modules

workflow METADATA {
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

process TRIMGALORE {

    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/trimmed", mode: 'copy', overwrite: true

    container 'quay.io/biocontainers/trim-galore:0.6.9--hdfd78af_0'

    input:
        tuple val(sample), path(fastq1), path(fastq2), path(gtf), path(fasta)

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
        --fastqc \
        ${fastq1} \
        ${fastq2}
    """
}

process FLASH {

    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/merged", mode: 'copy', overwrite: true

    container 'staphb/flash:1.2.11'

    input:
        tuple val(sample), path(fastq1), path(fastq2)

    output:
        tuple val(sample), path("${sample}.extendedFrags.fastq.gz"), path("${sample}.notCombined_1.fastq.gz"), path("${sample}.notCombined_2.fastq.gz"), emit: mergedfastq

    script:
    """
    flash \
      -m 18 \
      -x 0.25 \
      -t ${task.cpus} \
      --allow-outies \
      -z --compress-prog=gzip \
      --suffix=gz \
      -o ${sample} \
      ${fastq1} \
      ${fastq2}
    """
}

process FASTX {

    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/fastx", mode: 'copy', overwrite: true

    container 'quay.io/biocontainers/fastx_toolkit:0.0.14--hfc679d8_7'

    input:
        tuple val(sample), path(merge), path(unmerge1), path(unmerge2)
    
    output:
        tuple val(sample), path("${sample}.combined.fastq.gz"), emit: combinedfastq
        tuple val(sample), path("${sample}.combined.reverse.fastq.gz"), emit: reversefastq

    script:
    """
    # Reverse  complement unmerge2
    zcat ${unmerge2} | fastx_reverse_complement \
      -z \
      -o ${sample}.notCombined_2.fastq.revcomp.gz

    # Merge the files
    cat ${merge} ${unmerge1} ${sample}.notCombined_2.fastq.revcomp.gz > ${sample}.combined.fastq.gz

    # Reverse complement the combined file
    zcat ${sample}.combined.fastq.gz | fastx_reverse_complement \
    -z \
    -o ${sample}.combined.reverse.fastq.gz
    """
}

// Main pipeline

workflow {

    // Run the METADATA workflow
    METADATA(params.input)
    // METADATA.out.view { v -> "Channel is ${v}" }

    TRIMGALORE(METADATA.out)
    FLASH(TRIMGALORE.out.trimfastq)
    FASTX(FLASH.out.mergedfastq)

}


