#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process STAR_VIRAL_INDEX {

    tag "$fasta"
    label 'process_high'
    publishDir "${params.outdir}/indices", mode: 'copy', overwrite: true, pattern: '*.index'

    container 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_7'

    input:
        tuple path(gtf), path(fasta)
    
    output:
        tuple path(gtf), path(fasta), path("*.index"), emit: viralindex

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