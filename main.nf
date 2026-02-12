#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Include workflows and processes
include { METADATA } from './modules/shared/metadata.nf'
include { WORKFLOW_ABCD } from './workflows/workflow_abcd.nf'
include { WORKFLOW_POLYA } from './workflows/workflow_polya.nf'
include { STAR_VIRAL_INDEX } from './modules/shared/star_viral_index.nf'
include { STAR_JOINT_INDEX } from './modules/shared/star_joint_index.nf'
include { STAR_HOST } from './modules/shared/star_host.nf' 
include { SAMTOOLS_HOST } from './modules/shared/samtools_host.nf' 
include { STAR_VIRAL } from './modules/shared/star_viral.nf' 
include { SAMTOOLS_VIRAL } from './modules/shared/samtools_viral.nf' 
include { BEDTOOLS } from './modules/shared/bedtools.nf' 
include { QUANTIFICATION } from './modules/shared/quantification.nf'

// Main pipeline
workflow {
    
    // Run the METADATA workflow
    METADATA(params.input)

    // Run the pre-processing workflows
    WORKFLOW_ABCD(METADATA.out.branched_data.lib_abcd)
    WORKFLOW_POLYA(METADATA.out.branched_data.lib_polya)

    // Run STAR indexing    
    STAR_VIRAL_INDEX(METADATA.out.refs)
    STAR_JOINT_INDEX(METADATA.out.refs)

    // Key STAR_JOINT_INDEX outputs - to allow joining with pre-processing workflow outputs, ready for STAR_HOST
    joint_keyed = STAR_JOINT_INDEX.out.jointindex
        .map { gtf, fasta, jindex -> tuple("${gtf.getName().toString()}::${fasta.getName().toString()}", jindex) }
    //joint_keyed.view { "STAR jointindex keyed: ${it}" }

    // Mix pre-processing workflow outputs - for both libraries A,B,C,D and PolyA - into a single channel, and key by gtf+fasta to allow joining with STAR index channels
    all_preprocessed = WORKFLOW_ABCD.out.processed.mix(WORKFLOW_POLYA.out.processed) 
    //all_preprocessed.view { "All preprocessed: ${it}" }

    // Key all_preprocessed
    all_preprocessed_keyed = all_preprocessed
        .map { sample, combined, reverse, gtf, fasta, library -> tuple("${gtf.getName().toString()}::${fasta.getName().toString()}", sample, combined, reverse, fasta, gtf) }
    //all_preprocessed_keyed.view { "all_preprocessed keyed: ${it}" }

    // Combine channels and filter for matching keys
    joined_for_host = all_preprocessed_keyed.combine(joint_keyed, by: 0)  // Cartesian product of both channels, by first element (the key)
    // Optional: view to check
    //joined_for_host.view { "STAR_HOST input: ${it}" }

    // Run STAR_HOST
    STAR_HOST(joined_for_host)

    // Run SAMTOOLS_HOST
    SAMTOOLS_HOST(STAR_HOST.out.host_bam)

    // Key STAR_VIRAL_INDEX outputs
    viral_keyed = STAR_VIRAL_INDEX.out.viralindex
        .map { gtf, fasta, vindex -> tuple("${gtf.getName().toString()}::${fasta.getName().toString()}", vindex) }
    //viral_keyed.view { "STAR viralindex keyed: ${it}" }

    // Key SAMTOOLS_HOST viral outputs
    samtools_pre_keyed = SAMTOOLS_HOST.out.viral
        .map { sample, fasta, gtf, bam, fastq -> tuple("${gtf.getName().toString()}::${fasta.getName().toString()}", sample, fastq) }
    //samtools_pre_keyed.view { "SAMTOOLS_HOST viral keyed: ${it}" }

    // Combine channels and filter for matching keys
    joined_for_viral = samtools_pre_keyed.combine(viral_keyed, by: 0)  // Cartesian product of both channels, by first element (the key)
    //joined_for_viral.view { "STAR_VIRAL input: ${it}" }

    // Run STAR_VIRAL
    STAR_VIRAL(joined_for_viral)

    // Run SAMTOOLS_VIRAL
    SAMTOOLS_VIRAL(STAR_VIRAL.out.viral_bam)

    // Run BEDTOOLS
    BEDTOOLS(SAMTOOLS_VIRAL.out.bams)

    // Join SAMTOOLS_VIRAL outputs with METADATA.out.rawdata for downstream analysis
    joined_for_analysis = METADATA.out.rawdata.join(SAMTOOLS_VIRAL.out.bams)

    // Run QUANTIFICATION on the sorted viral BAMs
    QUANTIFICATION(
        joined_for_analysis,
        file("${projectDir}/bin/analysis.R"),
        file("${projectDir}/bin/functions.R")
    )

}