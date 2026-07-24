#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process FASTX {

    tag "$sample"
    label 'process_superhigh'

    container 'quay.io/biocontainers/fastx_toolkit:0.0.14--hfc679d8_7'

    input:
        tuple val(sample), path(fastq1), path(fastq2), path(gtf), path(fasta), val(library), path(merge), path(unmerge1), path(unmerge2)
    
    output:
        // Emit a single paired channel containing both combined and reverse fastq paths
        tuple val(sample), path("${sample}.combined.fastq.gz"), path("${sample}.combined.reverse.fastq.gz"), path(gtf), path(fasta), val(library), emit: fastx_pair

    script:
    """
    # Combine R1R2 collapsed reads with R1 singletones (LIB B or C or D), or R1R2 + R1 + R2 (LIB A)
    # R2 singletons are first reverse-complement to put them in the same orientation as R1 singletones.
    if test "${library}" == "A"
    then
        # Reverse complement unmerge2
        zcat ${unmerge2} | fastx_reverse_complement \
            -z \
            -o ${sample}.notCombined_2.fastq.revcomp.gz
        
        # Merge all 3 files (library A)
        cat ${merge} ${unmerge1} ${sample}.notCombined_2.fastq.revcomp.gz > ${sample}.combined.fastq.gz

        # Reverse complement the combined file
        zcat ${sample}.combined.fastq.gz | fastx_reverse_complement \
            -z \
            -o ${sample}.combined.reverse.fastq.gz
    else
        # For libraries B, C, D, we merge R1R2 with R1 singletons
        cat ${merge} ${unmerge1} > ${sample}.combined.fastq.gz

        # Reverse complement the combined file
        zcat ${sample}.combined.fastq.gz | fastx_reverse_complement \
            -z \
            -o ${sample}.combined.reverse.fastq.gz
    fi 
    """
    
}
