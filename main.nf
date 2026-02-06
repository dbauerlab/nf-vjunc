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

process STAR_HOST {

    tag "$sample"
    label 'process_high'
    publishDir "${params.outdir}/star_host", mode: 'copy', overwrite: true, pattern: '*.bam'

    container 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_7'

    input:
        tuple val(key), val(sample), path(combined), path(reverse), path(fasta), path(gtf), path(joint_index)

    output:
        tuple val(sample), path(fasta), path(gtf), path("${sample}_Aligned.out.bam"), emit: host_bam

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
        --outSAMtype BAM Unsorted \
        --outFileNamePrefix ${sample}_
    """

}

process SAMTOOLS_HOST {
    
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/samtools_host", mode: 'copy', overwrite: true, pattern: '*.bam'
    publishDir "${params.outdir}/samtools_host", mode: 'copy', overwrite: true, pattern: '*.fastq.gz'
    publishDir "${params.outdir}/samtools_host", mode: 'copy', overwrite: true, pattern: '*.idxstats'
    publishDir "${params.outdir}/samtools_host", mode: 'copy', overwrite: true, pattern: '*.flagstat'
    publishDir "${params.outdir}/samtools_host", mode: 'copy', overwrite: true, pattern: '*.coverage.txt'

    container 'quay.io/biocontainers/samtools:1.22.1--h96c455f_0'

    input:
        tuple val(sample), path(fasta), path(gtf), path(bam)

    output:
        tuple val(sample), path("${sample}.idxstats"), path("${sample}.flagstat"), path("${sample}.coverage.txt"), emit: stats
        tuple val(sample), path(fasta), path(gtf), path("${sample}.viral.bam"), path("${sample}.viral.fastq.gz"), emit: viral
        tuple val(sample), path("${sample}.unmapped.bam"), path("${sample}.unmapped.fastq.gz"), emit: unmapped

    script:
    """
    VIRAL_CHR=\$(grep '^>' ${fasta} | sed 's/ .*//' | sed 's/>//' | tr -d '[:space:]' | tr '\n' ' ' | tr '\r' ' ')

    # Sort, index, and generate stats on the pre-mapping BAM file
    samtools sort --threads $task.cpus -o ${sample}.sorted.bam $bam
    samtools index ${sample}.sorted.bam
    samtools idxstats ${sample}.sorted.bam > ${sample}.idxstats
    samtools flagstat ${sample}.sorted.bam > ${sample}.flagstat
    samtools depth -a -m 0 ${sample}.sorted.bam > ${sample}.coverage.txt

    # Collect reads mapping to viral genome into a new fastq file for further analysis
    samtools view --threads $task.cpus -b -o ${sample}.viral.bam ${sample}.sorted.bam \${VIRAL_CHR}
    samtools fastq --threads $task.cpus -0 ${sample}.viral.fastq.gz ${sample}.viral.bam

    # Collect unmapped reads into a new fastq file for further analysis
    samtools view --threads $task.cpus -b -f 4 -o ${sample}.unmapped.bam ${sample}.sorted.bam
    samtools fastq --threads $task.cpus -0 ${sample}.unmapped.fastq.gz ${sample}.unmapped.bam
    """
}

process STAR_VIRAL {

    tag "$sample"
    label 'process_high'
    publishDir "${params.outdir}/star_viral", mode: 'copy', overwrite: true, pattern: '*.bam'

    container 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_7'

    input:
        tuple val(key), val(sample), path(fastq), path(index)

    output:
        tuple val(sample), path("${sample}_Aligned.out.bam"), emit: viral_bam

    script:
    """
    STAR \
        --runThreadN $task.cpus \
        --genomeDir $index \
        --readFilesIn $fastq \
        --readFilesCommand zcat \
        --outFileNamePrefix ${sample}_ \
        --outReadsUnmapped Fastx \
        --outStd Log \
        --outSAMtype BAM Unsorted \
        --outSAMattributes Standard \
        --twopassMode Basic \
        --seedPerWindowNmax 30 \
        --alignIntronMin 1 \
        --outSJfilterOverhangMin 20 20 20 20 \
        --outSJfilterCountUniqueMin 1 1 1 1 \
        --outSJfilterCountTotalMin 1 1 1 1 \
        --outSJfilterDistToOtherSJmin 0 0 0 0 \
        --scoreGapNoncan 0 \
        --scoreGapGCAG 0 \
        --scoreGapATAC 0 \
        --alignSJoverhangMin 20 \
        --outFilterMatchNmin 40 \
        --outSJfilterReads All \
        --outSAMmultNmax 1 \
        --outFilterMismatchNoverLmax 0.1 \
        --alignEndsType Local \
        --outFilterType BySJout \
        --limitOutSJcollapsed 10000000 \
        --limitIObufferSize 1500000000 1500000000 \
        --alignSJstitchMismatchNmax 0 0 0 0 \
        --alignSJDBoverhangMin 20 \
        --alignSoftClipAtReferenceEnds Yes \
        --scoreGenomicLengthLog2scale 0 \
        --outFilterMultimapNmax 1
    """

}

process SAMTOOLS_VIRAL {
    
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/samtools_viral", mode: 'copy', overwrite: true, pattern: '*.bam'
    publishDir "${params.outdir}/samtools_viral", mode: 'copy', overwrite: true, pattern: '*.fastq.gz'
    publishDir "${params.outdir}/samtools_viral", mode: 'copy', overwrite: true, pattern: '*.idxstats'
    publishDir "${params.outdir}/samtools_viral", mode: 'copy', overwrite: true, pattern: '*.flagstat'
    publishDir "${params.outdir}/samtools_viral", mode: 'copy', overwrite: true, pattern: '*.txt'

    container 'quay.io/biocontainers/samtools:1.22.1--h96c455f_0'

    input:
        tuple val(sample), path(bam)

    output:
        tuple val(sample), path("${sample}.idxstats"), path("${sample}.flagstat"), path("${sample}.coverage.txt"), emit: stats
        tuple val(sample), path("${sample}.sorted.bam"), path("${sample}.spliced.bam"), emit: bams

    script:
    """
    # Sort, index, and generate stats on the pre-mapping BAM file
    samtools sort --threads $task.cpus -o ${sample}.sorted.bam $bam
    samtools index ${sample}.sorted.bam
    samtools idxstats ${sample}.sorted.bam > ${sample}.idxstats
    samtools flagstat ${sample}.sorted.bam > ${sample}.flagstat
    samtools depth -a -m 0 ${sample}.sorted.bam > ${sample}.coverage.txt

    # Collect reads that are spliced
    samtools view -h ${sample}.sorted.bam | awk -v OFS="\t" '\$0 ~ /^@/{print \$0;next;} \$6 ~ /N/' | samtools view -b -o ${sample}.spliced.bam
    samtools index ${sample}.spliced.bam
    """
}

process BEDTOOLS {
    
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/bedtools", mode: 'copy', overwrite: true, pattern: '*.sorted.bg'

    container 'quay.io/biocontainers/bedtools:2.31.1--h13024bc_3'

    input:
        tuple val(sample), path(sorted), path(spliced)

    output:
        tuple val(sample), path("${sample}.unstranded.cov.sorted.bg"), path("${sample}.plus.cov.sorted.bg"), path("${sample}.minus.cov.sorted.bg"), path("${sample}.spliced.unstranded.cov.sorted.bg"), path("${sample}.spliced.plus.cov.sorted.bg"), path("${sample}.spliced.minus.cov.sorted.bg"), emit: bedgraph
        
    script:
    """
    # Do coverage for reads from unfiltered BAM
    bedtools genomecov \
        -bga \
        -ibam ${sorted} > ${sample}.unstranded.cov.bg
    bedtools sort -i ${sample}.unstranded.cov.bg > ${sample}.unstranded.cov.sorted.bg
  
    bedtools genomecov \
        -bga \
        -strand "+" \
        -ibam ${sorted} > ${sample}.plus.cov.bg
    bedtools sort -i ${sample}.plus.cov.bg > ${sample}.plus.cov.sorted.bg

    bedtools genomecov \
        -bga \
        -strand "-" \
        -ibam ${sorted} > ${sample}.minus.cov.bg
    bedtools sort -i ${sample}.minus.cov.bg > ${sample}.minus.cov.sorted.bg
  
    # Do coverage for reads from spliced BAM
    bedtools genomecov \
        -bga \
        -ibam ${spliced} > ${sample}.spliced.unstranded.cov.bg
    bedtools sort -i ${sample}.spliced.unstranded.cov.bg > ${sample}.spliced.unstranded.cov.sorted.bg
  
    bedtools genomecov \
        -bga \
        -strand "+" \
        -ibam ${spliced} > ${sample}.spliced.plus.cov.bg
    bedtools sort -i ${sample}.spliced.plus.cov.bg > ${sample}.spliced.plus.cov.sorted.bg

    bedtools genomecov \
        -bga \
        -strand "-" \
        -ibam ${spliced} > ${sample}.spliced.minus.cov.bg
    bedtools sort -i ${sample}.spliced.minus.cov.bg > ${sample}.spliced.minus.cov.sorted.bg
    """
}

process R_ANALYSIS {

    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/r_analysis", mode: 'copy', overwrite: true
    container '/nemo/stp/babs/working/bootj/singularity/r-analysis.sif'

    input:
        tuple val(sample), path(fastq1), path(fastq2), path(gtf), path(fasta), val(library), path("${sample}.sorted.bam"), path("${sample}.spliced.bam")
        path analysis_script
        path functions_script
    
    script:
    """
    Rscript ${analysis_script} ${sample} ${sample}.sorted.bam ${gtf} ${fasta}
    """
}

// Main pipeline
workflow {
    
    // Run the METADATA workflow
    METADATA(params.input)

    // Run the pre-processing processes
    TRIMGALORE(METADATA.out.data)
    STAR_VIRAL_INDEX(METADATA.out.refs)
    STAR_JOINT_INDEX(METADATA.out.refs)
    UMITOOLS(TRIMGALORE.out.trimfastq)
    joined_for_hardtrim = METADATA.out.data.join(UMITOOLS.out.umifastq)
    HARDTRIM(joined_for_hardtrim)
    FLASH(HARDTRIM.out.clippedfastq)
    joined_for_fastx = METADATA.out.data.join(FLASH.out.mergedfastq)
    FASTX(joined_for_fastx)

    // Key STAR_JOINT_INDEX outputs
    joint_keyed = STAR_JOINT_INDEX.out.jointindex
        .map { gtf, fasta, jindex -> tuple("${gtf.getName().toString()}::${fasta.getName().toString()}", jindex) }
    //joint_keyed.view { "STAR jointindex keyed: ${it}" }

    // Key FASTX outputs
    fastx_keyed = FASTX.out.fastx_pair
        .map { sample, combined, reverse, gtf, fasta, library -> tuple("${gtf.getName().toString()}::${fasta.getName().toString()}", sample, combined, reverse, fasta, gtf) }
    //fastx_keyed.view { "FASTX fastx_pair keyed: ${it}" }

    // Combine channels and filter for matching keys
    joined_for_host = fastx_keyed.combine(joint_keyed, by: 0)  // Cartesian product of both channels, by first element (the key)
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
    // Optional: view to check
    //joined_for_viral.view { "STAR_VIRAL input: ${it}" }

    // Run STAR_VIRAL
    STAR_VIRAL(joined_for_viral)

    // Run SAMTOOLS_VIRAL
    SAMTOOLS_VIRAL(STAR_VIRAL.out.viral_bam)

    // Run BEDTOOLS
    BEDTOOLS(SAMTOOLS_VIRAL.out.bams)

    // Join SAMTOOLS_VIRAL outputs with METADATA.out.data for downstream analysis
    joined_for_analysis = METADATA.out.data.join(SAMTOOLS_VIRAL.out.bams)

    // Run R_ANALYSIS on the sorted viral BAMs
    R_ANALYSIS(
        joined_for_analysis,
        file("${projectDir}/R/analysis.R"),
        file("${projectDir}/R/functions.R")
    )

}