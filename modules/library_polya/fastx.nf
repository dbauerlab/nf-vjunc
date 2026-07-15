#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process FASTX {

    tag "$sample"
    label 'process_medium'

    container 'quay.io/biocontainers/fastx_toolkit:0.0.14--hfc679d8_7'

    input:
        tuple val(sample), path(fastq1), path(fastq2), path(gtf), path(fasta), val(library), path(merge), path(unmerge1), path(unmerge2)
    
    output:
        // Emit a single paired channel containing both combined and reverse fastq paths
        tuple val(sample), path("${sample}.combined.fastq.gz"), path("${sample}.combined.reverse.fastq.gz"), path(gtf), path(fasta), val(library), emit: fastx_pair

    script:
    """
    # Merge R1+R2 collapsed reads with R1 singletons (combined read is on R1 strand so no need to edit R1 singletons)
    # Also merge R1+R2 collapsed reads with R2 singletons (R2 singletons need to be reverse complemented)

    # Reverse complement unmerge2
    zcat ${unmerge2} | fastx_reverse_complement \
        -z \
        -o ${sample}.notCombined_2.fastq.revcomp.gz
        
    # Merge all 3 files
    cat ${merge} ${unmerge1} ${sample}.notCombined_2.fastq.revcomp.gz > ${sample}.combined.fastq.gz

    # Reverse complement - everything is now on R1 strand (which is negative strand), so reverse complement to get positive strand
    zcat ${sample}.combined.fastq.gz | fastx_reverse_complement \
        -z \
        -o ${sample}.combined.reverse.fastq.gz
    """
    
}