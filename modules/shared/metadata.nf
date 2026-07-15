#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Modules

workflow METADATA {
    take: csv
    main:
        // Parse CSV and group rows by sample name, collecting fastq files across lanes into lists
        Channel
            .fromPath( csv )
            .splitCsv(header:true)
            .map { row -> [ row.sample, 
                            file(row.fastq1, checkIfExists: true),
                            file(row.fastq2, checkIfExists: true),
                            file(row.gtf, checkIfExists: true),
                            file(row.fasta, checkIfExists: true),
                            row.library ]  }
            .groupTuple(by: 0)
            .map { sample, fastq1s, fastq2s, gtfs, fastas, libraries ->
                // gtf, fasta, library must be identical across lanes for the same sample
                [ sample, fastq1s, fastq2s, gtfs[0], fastas[0], libraries[0] ]
            }
            .set { grouped_data }

        // rawdata: single representative fastq per sample (first lane) for downstream metadata joins
        grouped_data
            .map { sample, fastq1s, fastq2s, gtf, fasta, library ->
                [ sample, fastq1s[0], fastq2s[0], gtf, fasta, library ]
            }
            .set { rawdata }

        // Branch grouped_data (fastq1s/fastq2s are lists) by library type for pre-processing workflows
        grouped_data
            .branch {
                lib_abcd: it[5] in ['A', 'B', 'C', 'D']
                    return it // Returns tuple: (sample_id, [fastq1s], [fastq2s], gtf, fasta, library)
                lib_polya: it[5] == 'PolyA' 
                    return it // Returns tuple: (sample_id, [fastq1s], [fastq2s], gtf, fasta, library)
            } 
            .set { branched_data }

        // create refs channel with unique fasta+gtf pairs
        grouped_data
            .map { sample, fastq1s, fastq2s, gtf, fasta, library ->
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