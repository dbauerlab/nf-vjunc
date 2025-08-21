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
            .set { data }
    emit:
        data
}

process TRIMGALORE {

    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/adapter_trim", mode: 'copy', overwrite: true

    container 'quay.io/biocontainers/trim-galore:0.6.9--hdfd78af_0'

    input:
        tuple val(sample), path(fastq1), path(fastq2), path(gtf), path(fasta), val(library)

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
        ${fastq1} \
        ${fastq2}
    """
}

process UMITOOLS {

    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/umitools", mode: 'copy', overwrite: true

    container 'quay.io/biocontainers/umi_tools:1.1.5--py39hbcbf7aa_4'

    input:
        tuple val(sample), path(fastq1), path(fastq2)

    output:
        tuple val(sample), path("${sample}.umi_extracted.R1.fastq.gz"), path("${sample}.umi_extracted.R2.fastq.gz"), emit: umifastq
        tuple val(sample), path("${sample}.umi_log"), emit: umilog

    script:
    """
    umi_tools extract \
        -I ${sample}_val_1.fq.gz \
        --read2-in=${sample}_val_2.fq.gz \
        --bc-pattern="NNNNNNNNNN" \
        --bc-pattern2="NNNNNNNNNN" \
        --stdout=${sample}.umi_extracted.R1.fastq.gz \
        --read2-out=${sample}.umi_extracted.R2.fastq.gz \
        --log=${sample}.umi_log
    """
}

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
        tuple val(sample), path(fastq1), path(fastq2), path(gtf), path(fasta), val(library), path(merge), path(unmerge1), path(unmerge2)
    
    output:
        tuple val(sample), path("${sample}.combined.fastq.gz"), emit: combinedfastq
        tuple val(sample), path("${sample}.combined.reverse.fastq.gz"), emit: reversefastq

    script:
    """
    # Combine R1R2 collapsed reads with R1 singletones (LIB B or C or D), or R1R2 + R1 + R2 (LIB A)
    # R2 singletons are first reverse-complement to put them in the same orientation as R1 singletones.
    if test "${LIB}" == "A"
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

// Main pipeline

workflow {

    // Run the METADATA workflow
    METADATA(params.input)
    // METADATA.out.view { v -> "Channel is ${v}" }

    TRIMGALORE(METADATA.out)
    UMITOOLS(TRIMGALORE.out.trimfastq)
    joined_for_hardtrim = METADATA.out.join(UMITOOLS.out.umifastq)
    HARDTRIM(joined_for_hardtrim)
    FLASH(HARDTRIM.out.clippedfastq)
    joined_for_fastx = METADATA.out.join(FLASH.out.mergedfastq)
    FASTX(joined_for_fastx)

}


