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

        // create refs channel with unique fasta+gtf pairs
        data
            .map { sample, fastq1, fastq2, gtf, fasta, library ->
                tuple(gtf, fasta)
            }
            .unique()
            .set { refs }
    emit:
        data
        refs
}

process STAR_VIRAL_INDEX {

    tag "$fasta"
    label 'process_high'
    publishDir "${params.outdir}/indices", mode: 'copy', overwrite: true, pattern: '*.index'

    container 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_7'

    input:
        tuple path(gtf), path(fasta)
    
    output:
        tuple path(fasta), path(gtf), path("*.index"), emit: viralindex

    script:
    """
    # Compute a safe prefix from the fasta filename
    fname=\$(basename "$fasta")
    prefix=\${fname%%.*}   # strips first extension, handles e.g. hg38.fa or hg38.fa.gz

    mkdir -p "\${prefix}.index"

    # compute genome length (sum of sequence lengths). Use awk to avoid extra tools.
    genlen=\$(awk '/^>/ { if (seqlen){ total+=seqlen; seqlen=0 } ; next } { seqlen+=length(\$0) } END { total+=seqlen; print total }' $fasta)

    # compute genomeSAindexNbases using recommended heuristic:
    # genomeSAindexNbases = min(14, max(4, int(log2(genomeLength)/2 - 1)))
    saIndex=\$(awk -v L="\${genlen}" 'BEGIN{ if(L<=0){print 14; exit} g=log(L)/log(2); v=int(g/2 - 1); if(v<4) v=4; if(v>14) v=14; print v }')

    echo "Genome length: \${genlen}, using genomeSAindexNbases=\${saIndex}"

    STAR \
        --runMode genomeGenerate \
    --genomeDir \${prefix}.index/ \
        --genomeFastaFiles $fasta \
        --sjdbGTFfile $gtf \
        --runThreadN $task.cpus \
    --genomeSAindexNbases \${saIndex}
    """
}

process STAR_JOINT_INDEX {

    tag "$fasta"
    label 'process_high'
    publishDir "${params.outdir}/indices", mode: 'copy', overwrite: true, pattern: '*.index'
    publishDir "${params.outdir}/joint_fasta", mode: 'copy',  overwrite: true, pattern: '*.fa'
    publishDir "${params.outdir}/joint_gtf", mode: 'copy',  overwrite: true, pattern: '*.gtf'

    container 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_7'

    input:
        tuple path(gtf), path(fasta)
    
    output:
        tuple path(gtf), path(fasta), path("*combined.fa"), path("*combined.gtf"), path("*.index"), emit: jointindex

    script:
    """
    # Validate host files are provided
    if test -z "${params.host_fasta}" || test -z "${params.host_gtf}" ; then
        echo "ERROR: params.host_fasta and params.host_gtf must be set for STAR_JOINT_INDEX" >&2
        exit 2
    fi

    # Compute a safe prefix from the sample fasta file name
    sfname=\$(basename "$fasta")
    sprefix=\${sfname%%.*}

    # derive host basename
    hfname=\$(basename "${params.host_fasta}")
    hprefix=\${hfname%%.*}

    joint_prefix=\${sprefix}_\${hprefix}
    mkdir -p "\${joint_prefix}.joint.index"

    # create combined fasta and gtf
    combined_fasta=\${joint_prefix}.combined.fa
    combined_gtf=\${joint_prefix}.combined.gtf

    # helper: cat or gunzip -c depending on file extension
    cat_or_zcat() {
        f="\$1"
        case "\$f" in
            *.gz) gunzip -c "\$f" ;; 
            *) cat "\$f" ;;
        esac
    }

    # concatenate host then sample (order: host then sample)
    cat_or_zcat "${params.host_fasta}" > \${combined_fasta}
    cat_or_zcat "$fasta" >> \${combined_fasta}

    cat_or_zcat "${params.host_gtf}" > \${combined_gtf}
    cat_or_zcat "$gtf" >> \${combined_gtf}

    # compute genome length (sum of sequence lengths). Use awk to avoid extra tools.
    genlen=\$(awk '/^>/ { if (seqlen){ total+=seqlen; seqlen=0 } ; next } { seqlen+=length(\$0) } END { total+=seqlen; print total }' \${combined_fasta})

    # compute genomeSAindexNbases using recommended heuristic:
    saIndex=\$(awk -v L="\${genlen}" 'BEGIN{ if(L<=0){print 14; exit} g=log(L)/log(2); v=int(g/2 - 1); if(v<4) v=4; if(v>14) v=14; print v }')

    echo "Joint genome length: \${genlen}, using genomeSAindexNbases=\${saIndex}"

    STAR \
        --runMode genomeGenerate \
        --genomeDir \${joint_prefix}.joint.index/ \
        --genomeFastaFiles \${combined_fasta} \
        --sjdbGTFfile \${combined_gtf} \
        --runThreadN $task.cpus \
        --genomeSAindexNbases \${saIndex}
    """
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

process STAR_PREMAP {

    tag "$sample"
    label 'process_high'
    publishDir "${params.outdir}/star_premap", mode: 'copy', overwrite: true

    container 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_7'

    input:
        tuple val(sample), path(combined), path(reverse), path(gtf), path(fasta), val(library), path(joint_combined_fa), path(joint_combined_gtf), path(joint_index)

    output:
        tuple val(sample), path(bam), emit: premap_bam

    script:
    """
    STAR \
        --runThreadN $task.cpus \
        --genomeDir $joint_index \
        --readFilesIn $reverse \
        --readFilesCommand zcat \
        --twopassMode Basic \
        --outReadsUnmapped None \
        --outSAMunmapped Within \
        --outSAMtype BAM Unsorted
    """

}


// Main pipeline

workflow {

    // Run the METADATA workflow
    METADATA(params.input)
    // METADATA.out.view { v -> "Channel is ${v}" }

    TRIMGALORE(METADATA.out.data)
    STAR_VIRAL_INDEX(METADATA.out.refs)
    STAR_JOINT_INDEX(METADATA.out.refs)
    UMITOOLS(TRIMGALORE.out.trimfastq)
    joined_for_hardtrim = METADATA.out.data.join(UMITOOLS.out.umifastq)
    HARDTRIM(joined_for_hardtrim)
    FLASH(HARDTRIM.out.clippedfastq)
    joined_for_fastx = METADATA.out.data.join(FLASH.out.mergedfastq)
    FASTX(joined_for_fastx)

    // Build joint index keyed by original gtf::fasta so we can match it to samples
    joint_keyed = STAR_JOINT_INDEX.out.jointindex
        .map{ gtf, fasta, combined_fa, combined_gtf, idx -> tuple("${gtf.toString()}::${fasta.toString()}", gtf, fasta, combined_fa, combined_gtf, idx) }

    // Diagnostic: show joint_keyed contents if available
    if (binding.hasVariable('joint_keyed') && joint_keyed) {
        joint_keyed.view { v -> "joint_keyed: ${v}" }
    } else {
        println "Warning: joint_keyed channel is null/empty"
    }

    // Pair FASTX combined and reverse outputs per sample
    fastx_pairs = FASTX.out.combinedfastq
        .join(FASTX.out.reversefastq)
        .map{ left, right ->
            sample = left[0]
            combined = left[1]
            reverse = right[1]
            tuple(sample, combined, reverse)
        }

    // Diagnostic: show fastx_pairs contents if available
    if (binding.hasVariable('fastx_pairs') && fastx_pairs) {
        fastx_pairs.view { v -> "fastx_pairs: ${v}" }
    } else {
        println "Warning: fastx_pairs channel is null/empty"
    }

    // Attach sample metadata (gtf,fasta,library) to FASTX outputs and key by gtf::fasta
    fastx_with_meta = fastx_pairs
        .join( METADATA.out.data.map{ sample, fastq1, fastq2, gtf, fasta, library -> tuple(sample, gtf, fasta, library) } )
        .map{ left, right ->
            sample = left[0]
            combined = left[1]
            reverse = left[2]
            gtf = right[1]
            fasta = right[2]
            library = right[3]
            key = "${gtf.toString()}::${fasta.toString()}"
            tuple(key, sample, combined, reverse, gtf, fasta, library)
        }

    // Diagnostic: show fastx_with_meta contents if available
    if (binding.hasVariable('fastx_with_meta') && fastx_with_meta) {
        fastx_with_meta.view { v -> "fastx_with_meta: ${v}" }
    } else {
        println "Warning: fastx_with_meta channel is null/empty"
    }

    // Join FASTX+meta with joint index so each sample receives the matching joint index path
    aligned_inputs = fastx_with_meta
        .join(joint_keyed)
        .map{ left, right ->
            // left: [key, sample, combined, reverse, gtf, fasta, library]
            // right: [key, gtf, fasta, combined_fa, combined_gtf, idx]
            sample = left[1]
            combined = left[2]
            reverse = left[3]
            gtf = left[4]
            fasta = left[5]
            library = left[6]
            joint_combined_fa = right[3]
            joint_combined_gtf = right[4]
            joint_index = right[5]
            tuple(sample, combined, reverse, gtf, fasta, library, joint_combined_fa, joint_combined_gtf, joint_index)
        }
        .set { aligned_inputs }
    
    // Print a single collected summary after upstream closes. collect() buffers all items in memory.
    if (aligned_inputs) {
        aligned_inputs
            .collect()
            .view { v -> "Aligned inputs (all): ${v}" }
    } else {
        println "Warning: aligned_inputs channel is null/empty; skipping collect/view"
    }

    }
