#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Include processes needed for this workflow
include { CAT_FASTQ } from '../modules/shared/cat_fastq.nf'
include { TRIMGALORE } from '../modules/library_abcd/trimgalore.nf'
include { UMITOOLS } from '../modules/library_abcd/umitools.nf'
include { HARDTRIM } from '../modules/library_abcd/hardtrim.nf'
include { FLASH } from '../modules/library_abcd/flash.nf'
include { FASTX } from '../modules/library_abcd/fastx.nf'

workflow WORKFLOW_ABCD {
    take:
        data  // tuple: sample_id, fastq1, fastq2, gtf, fasta, library
    
    main:
        // Concatenate across lanes first; single-lane samples are copied as-is
        CAT_FASTQ(data)
        // Run the pre-processing processes
        TRIMGALORE(CAT_FASTQ.out.catfastq)
        UMITOOLS(TRIMGALORE.out.trimfastq)
        joined_for_hardtrim = CAT_FASTQ.out.catfastq.join(UMITOOLS.out.umifastq)
        HARDTRIM(joined_for_hardtrim)
        FLASH(HARDTRIM.out.clippedfastq)
        joined_for_fastx = CAT_FASTQ.out.catfastq.join(FLASH.out.mergedfastq)
        results_ch = FASTX(joined_for_fastx)
    
    emit:
        processed = results_ch
    
}