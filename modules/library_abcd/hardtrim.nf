#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process HARDTRIM {

    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/hardtrim", mode: 'copy', overwrite: true

    container 'quay.io/biocontainers/fastx_toolkit:0.0.14--hfc679d8_7'

    input:
        tuple val(sample), path(fastq1), path(fastq2), path(gtf), path(fasta), val(library), path(umifastq1), path(umifastq2)

    output:
        tuple val(sample), path("${sample}.collapseReady.R1.fastq.gz"), path("${sample}.collapseReady.R2.fastq.gz"), emit: clippedfastq

    script:
    """
    ## Hard-clip the PCR primer from the sequences downstream of the removed R1 UMI. The amount of the sequence to clip is library specific. For library B the upper limit of the expected size is 27bp.

    if test "${library}" != "A"
    then 
        HARDCLIP=0
        if test "${library}" = "B"
        then
            HARDCLIP=27
        fi
        if test "${library}" = "C"
        then
            HARDCLIP=19
        fi
        if test "${library}" = "D"
        then
            HARDCLIP=19
        fi
        HARDCLIP=`expr \${HARDCLIP} + 1`
        # Run the hard-trim on R1
        zcat ${umifastq1} | fastx_trimmer \
            -z \
            -f "\${HARDCLIP}" \
            -o ${sample}.collapseReady.R1.fastq.gz
        # Copy R2 as it is
        cp ${umifastq2} ${sample}.collapseReady.R2.fastq.gz
    else
        # If library is A, just copy the files
        cp ${umifastq1} ${sample}.collapseReady.R1.fastq.gz
        cp ${umifastq2} ${sample}.collapseReady.R2.fastq.gz
    fi
    """
    
}