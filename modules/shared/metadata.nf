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
                            file(row.fasta, checkIfExists: true),
                            row.library ]  }
            .set { rawdata }

        // branch rawdata depending on library type - if A,B,C,D, send to WORKFLOW_ABCD, if PolyA, send to WORKFLOW_POLYA 
        rawdata
            .branch {
                lib_abcd: it[5] in ['A', 'B', 'C', 'D']
                    return it // Returns full tuple: (sample_id, fastq1, fastq2, gtf, fasta, library)
                lib_polya: it[5] == 'PolyA' 
                    return it // Returns full tuple: (sample_id, fastq1, fastq2, gtf, fasta, library)
            } 
            .set { branched_data }

        // create refs channel with unique fasta+gtf pairs
        rawdata
            .map { sample, fastq1, fastq2, gtf, fasta, library ->
                tuple(gtf, fasta)
            }
            .unique()
            .set { refs }
    emit:
        rawdata
        lib_abcd = branched_data.lib_abcd
        lib_polya = branched_data.lib_polya
        refs

}