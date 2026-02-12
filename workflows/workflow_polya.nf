#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Include processes needed for this workflow
include { TRIMGALORE } from '../modules/library_polya/trimgalore.nf'
include { FLASH } from '../modules/library_polya/flash.nf'
include { FASTX } from '../modules/library_polya/fastx.nf'

workflow WORKFLOW_POLYA {
    take:
        data  // tuple: sample_id, fastq1, fastq2, gtf, fasta, library
    
    main:
        // Run the pre-processing processes
        TRIMGALORE(data)
        FLASH(TRIMGALORE.out.trimfastq)
        joined_for_fastx = data.join(FLASH.out.mergedfastq)
        results_ch = FASTX(joined_for_fastx)
    
    emit:
        processed = results_ch
    
}