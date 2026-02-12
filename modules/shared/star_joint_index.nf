#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process STAR_JOINT_INDEX {

    tag "$fasta"
    label 'process_high'
    publishDir "${params.outdir}/indices", mode: 'copy', overwrite: true, pattern: '*.index'

    container 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_7'

    input:
        tuple path(gtf), path(fasta)
    
    output:
        tuple path(gtf), path(fasta), path("*.index"), emit: jointindex

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