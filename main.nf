#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

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

workflow {

    // Run the METADATA workflow
    METADATA(params.input)
    METADATA.out.view { v -> "Channel is ${v}" }

}


