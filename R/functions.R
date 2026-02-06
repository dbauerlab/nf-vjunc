# Functions for RNA recombination/junction/splicing analysis
# The functions have been altered from the original versions to work on a by sample basis, 
# taking the aligned sorted BAM file and the sample matches GTF as input and then outputting tables of junctions and segments for each sample which can then be combined for downstream analysis. 
# The sample BAM and GTF will be bound together in nextflow and then passed together to this script for processing.

################################################################################
#                                 Process BAM                                  #
################################################################################

ProcessBAMs <- function(sample, bam, gtf) { 
    # Define output list
    results.list <- list()
    # Output dataframe and columns
    xtab <- data.frame(
      total_reads = NA,
      spliced_reads = NA,
      positive_strand = NA,
      negative_strand = NA,
      genome_reads = NA
    )
    # Read in all reads from BAM file
    message('Reading BAM...')
    single.ga <- readGAlignments(
      file = BamFile(bam),
      use.names = TRUE,
      param = ScanBamParam(
        what = c("seq"),
        flag = scanBamFlag(
          isSecondMateRead = NA
        )
      ),
      with.which_label = FALSE
    )
    # Count total number of reads
    message('Counting total number of reads...')
    xtab$total_reads[1] <- length(single.ga)
    # Count number of +ve and -ve strand reads
    message('Counting +ve & -ve reads...')
    strand_info <- strand(single.ga)
    xtab$positive_strand[1] <- sum(strand_info == "+")
    xtab$negative_strand[1] <- sum(strand_info == "-")     
    # Limit by strand
    single.ga <- single.ga[strand(single.ga) %in% "+", ]  
    # Isolate spliced and unspliced
    message('Splitting spliced and unspliced...')
    spliced.ga <- single.ga[grepl("N", cigar(single.ga)), ]
    unspliced.ga <- single.ga[!grepl("N", cigar(single.ga)), ]  
    # Count number of spliced and unspliced reads
    message('Counting spliced and unspliced...')
    xtab$spliced_reads[1] <- length(spliced.ga)
    xtab$unspliced_reads[1] <- length(unspliced.ga) 
    # Count number of genome reads (for background)
    message('Counting genome reads...')
    # Define the genome background region as +/-500bp from the start of the S sgRNA intron
    S.Start <- max(gtf[gtf$gene_id == 'S']@ranges@start)
    ORF1ab.Start <- gtf[gtf$gene_id == 'ORF1ab' & gtf$type == 'exon']@ranges@start
    midpoint <- round((S.Start + ORF1ab.Start) / 2)
    genomeRegion <- GRanges(gtf@seqnames@values,
                            IRanges(start = midpoint, width = 1000),
                            strand = "+")
    # Count unspliced reads that overlap
    genome_bg.ga <- single.ga[single.ga %over% genomeRegion & !grepl("N", cigar(single.ga)), ]
    xtab$genome_reads[1] <- length(genome_bg.ga) 
    # Segment counts
    message('Calculating segment counts...')
    total.spliced.cnt <- table(seqnames(single.ga))[seqlevels(gtf)]
    seg.spliced.cnt <- table(seqnames(spliced.ga))[seqlevels(gtf)]
    seg.df <- data.frame(
      segment = seqlevels(gtf),
      total   = as.numeric(total.spliced.cnt),
      spliced = as.numeric(seg.spliced.cnt)
    )
    # Summarise junctions
    sample_junc.gr <- summarizeJunctions(single.ga, with.revmap = FALSE)
    # Add results to list
    results.list$xtab <- xtab
    results.list$segment_cnt <- seg.df
    results.list$junctions <- sample_junc.gr
    # Save xtab for sample as csv
    write.csv(results.list$xtab, file = paste0(sample, '_stats.csv'))
    # Save segment counts for sample as csv
    write.csv(results.list$segment_cnt, file = paste0(sample, '_segment_counts.csv'))
    # Save junctions for sample as RDS
    saveRDS(results.list$junctions, file = paste0(sample, '_junction_obj.RDS'))
    # Return results
    return(results.list)
}

################################################################################
#                           Table of all junctions                             #
################################################################################

# Compile a table of all junctions with their depth of coverage.
JunctionTable <- function(sample, results.list) {
    # Empty list for downstream
    sample_junc.list <- list()
    # Get junctions for sample
    sample_junc.gr <- results.list$junctions
    # If no junctions, create empty dataframe with 0 depth
    if (length(sample_junc.gr) == 0) {
      sample_junc.df <- data.frame(
        Sample_ID         = sample,
        chr            = NA,
        donor_site     = -1,
        acceptor_site  = -1,
        junction_depth = 0,
        stringsAsFactors = F
      )
    } else {
      # Create dataframe of junctions with depth
      sample_junc.df <- data.frame(
        Sample_ID      = sample,
        chr            = seqnames(sample_junc.gr),
        donor_site     = start(sample_junc.gr),
        acceptor_site  = end(sample_junc.gr),
        junction_depth = sample_junc.gr$score,
        stringsAsFactors = F
      )
    }
    # Junction frequency = junction depth / total junction depth
    sample_junc.df$total_junction_depth = sum(sample_junc.df$junction_depth)
    sample_junc.df$junction_freq        = sample_junc.df$junction_depth / sample_junc.df$total_junction_depth
    sample_junc.df$log10_junction_freq  = log10(sample_junc.df$junction_freq)
    # Replace NA and NaN with 0
    sample_junc.df$total_junction_depth[is.na(sample_junc.df$total_junction_depth) |
                                          is.nan(sample_junc.df$total_junction_depth)] <- 0
    sample_junc.df$junction_freq[is.na(sample_junc.df$junction_freq) |
                                   is.nan(sample_junc.df$junction_freq)] <- 0
    sample_junc.df$log10_junction_freq[is.na(sample_junc.df$log10_junction_freq) |
                                         is.nan(sample_junc.df$log10_junction_freq)] <- 0
    # Add gap size
    # Remove junction depth = 0
    sample_junc.df <- sample_junc.df %>%
      mutate(gap_size = acceptor_site - donor_site) %>%
      filter(junction_depth > 0)
    # Save
    write.csv(sample_junc.df, file = paste0(sample, '_junction_table.csv'))
    return(sample_junc.df)
}

################################################################################
#                      Add sequence of donor/acceptor                          #
################################################################################

AddJuncSeq <- function(sample, junc.tab, fasta){
  # Check if junction table is empty or has no valid junctions
  if (nrow(junc.tab) == 0 || all(is.na(junc.tab$chr)) || all(junc.tab$junction_depth == 0)) {
    message("No valid junctions found for sample ", sample, ". Skipping sequence addition.")
    # Return empty or minimal dataframe
    junc.seq.dat <- junc.tab
    junc.seq.dat$donor_seq <- character(0)
    junc.seq.dat$acceptor_seq <- character(0)
    # Save empty results
    saveRDS(junc.seq.dat, file = paste0(sample, '_junction_seq.RDS'))
    write.csv(junc.seq.dat, file = paste0(sample, '_junction_seq.csv'))
    return(junc.seq.dat)
  }
  # Index fasta
  indexFa(fasta)
  # Load FASTA file
  fastaObj <- FaFile(fasta)
  # Function to get donor/acceptor sequences for a given subset of junc.dat
  get_junc_seqs <- function(subdat, fa) {
    # Donor GRanges
    donor.gr <- GRanges(
      seqnames = subdat$chr,
      ranges   = IRanges(start = subdat$donor_site - 10, 
                         end   = subdat$donor_site + 10),
      strand = "+"
    )
    # Acceptor GRanges
    acceptor.gr <- GRanges(
      seqnames = subdat$chr,
      ranges   = IRanges(start = subdat$acceptor_site - 10, 
                         end   = subdat$acceptor_site + 10),
      strand = "+"
    )
    donor.seqs <- getSeq(fa, donor.gr)
    acceptor.seqs <- getSeq(fa, acceptor.gr)
    subdat$donor_seq    <- as.character(donor.seqs)
    subdat$acceptor_seq <- as.character(acceptor.seqs)
    subdat
  }
  # Get seqs
  junc.seq.dat <- get_junc_seqs(junc.tab, fastaObj)
  # Sort by junction freq
  junc.seq.dat <- junc.seq.dat %>%
    arrange(desc(log10_junction_freq))
  # Save for sample
  saveRDS(junc.seq.dat, file = paste0(sample, '_junction_seq.RDS'))
  write.csv(junc.seq.dat, file = paste0(sample, '_junction_seq.csv'))
  # Return
  return(junc.seq.dat)
}

################################################################################
#                             Junction Class                                   #
################################################################################

## Define classes of spliced junction reads based on overlap with canonical mRNA.
## The junction reads are redefined as the "gaps" (introns) between the exons, to match the gtf definitions.
## canonical = matches 2/2 junctions
## alternative = matches 1/2 junctions
## defective = matches 0/2 junctions

JuncClass <- function(sample, junc.seq.dat, gtf) {
    # Define canonical junctions from the gtf file(s).
    canonical_junctions.gr <- list()
    gtf.gr <- gtf
    # Get exons for each transcript and then get the gaps (introns) between them as the canonical junctions
    canonical.gr <- lapply(split(gtf.gr, gtf.gr$gene_id), function(x) {
      x[x$type %in% "exon", ]
    })
    canonical.gr <- canonical.gr[sapply(canonical.gr, length) > 1]      # Only splicing sgRNA
    # Loop through genes and transcripts to get the gaps (introns) as the canonical junctions
    for (g in 1:length(canonical.gr)) {
      canonical.tx.gr <- split(canonical.gr[[g]], canonical.gr[[g]]$transcript_id)
      canonical.tx.gr <- canonical.tx.gr[lapply(canonical.tx.gr, length) >
                                           1]
      # Loop through transcripts to get the gaps (introns) as the canonical junctions
      for (t in 1:length(canonical.tx.gr)) {
        canonical.tx.gaps.gr <- gaps(canonical.tx.gr[[t]])
        canonical.tx.gaps.gr$gene_id <- names(canonical.gr)[g]
        canonical.tx.gaps.gr$transcript_name <- names(canonical.tx.gr)[t]
        canonical_junctions.gr <- c(canonical_junctions.gr, canonical.tx.gaps.gr)
      }
    }
    # Unlist and remove duplicates
    canonical_junctions.gr <- unlist(GRangesList(canonical_junctions.gr))
    canonical_junctions.gr$key <- paste(
      seqnames(canonical_junctions.gr),
      start(canonical_junctions.gr),
      end(canonical_junctions.gr),
      sep = '.'
    )
    canonical_junctions.gr <- canonical_junctions.gr[!duplicated(canonical_junctions.gr$key), ]
    ## Classify junctions
    # Default is defective
    junc.seq.dat$class <- factor("defective",
                                 levels = c("canonical", "alternative", "defective"))
    junc.seq.dat$gene_id <- ""
    junc.seq.dat$transcript_name <- ""
    # Loop through table and classify as either canonical or alternative
    for (r in 1:nrow(junc.seq.dat)) {
      # Canonical
      canonical.idx <-
        seqnames(canonical_junctions.gr) %in% junc.seq.dat$chr[r] &
        (
          start(canonical_junctions.gr) %in% junc.seq.dat$donor_site[r] &
            end(canonical_junctions.gr) %in% junc.seq.dat$acceptor_site[r]
        )
      if (any(canonical.idx)) {
        junc.seq.dat$class[r] <- "canonical"
        junc.seq.dat$gene_id[r] <- canonical_junctions.gr[canonical.idx, ]$gene_id
        junc.seq.dat$transcript_name[r] <- canonical_junctions.gr[canonical.idx, ]$transcript_name
      } else {
        # Alternative
        alternative.idx <-
          seqnames(canonical_junctions.gr) %in% junc.seq.dat$chr[r] &
          (((
            start(canonical_junctions.gr) %in% junc.seq.dat$donor_site[r]
          ) &
            (
              !end(canonical_junctions.gr) %in% junc.seq.dat$acceptor_site[r]
            )
          ) |
            ((
              !start(canonical_junctions.gr) %in% junc.seq.dat$donor_site[r]
            ) &
              (
                end(canonical_junctions.gr) %in% junc.seq.dat$acceptor_site[r]
              )
            ))
        if (any(alternative.idx)) {
          junc.seq.dat$class[r] <- "alternative"
          junc.seq.dat$gene_id[r] <- paste(unique(canonical_junctions.gr[alternative.idx, ]$gene_id), collapse =
                                             ';')
          junc.seq.dat$transcript_name[r] <- paste(unique(canonical_junctions.gr[alternative.idx, ]$transcript_name),
                                                   collapse = ';')
        }
      }
    }
    # Save for sample
    saveRDS(junc.seq.dat, file = paste0(sample, '_junction_class.RDS'))
    write.csv(junc.seq.dat,
              file = paste0(sample, '_junction_class.csv'))
    # Return
    return(junc.seq.dat)
}