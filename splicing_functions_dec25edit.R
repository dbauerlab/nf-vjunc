# Functions for splicing analysis

################################################################################
#                                  Import GTF                                  #
################################################################################

ImportGTF <- function(dir){
  # Find gtf files
  gtfFiles <- list.files(path = dir,
                         pattern = "\\.gtf$",
                         recursive = T,
                         full.names = T)
  # Get genome names 
  genomes <- sub(".*\\/([^\\/]+)$", "\\1", gtfFiles)
  genomes <- gsub('.gtf', '', genomes)
  # Import gtfs
  gtf.dat.list <- lapply(gtfFiles, function(x){
    rtracklayer::import(x)
  })
  # Define names of gtf list
  names(gtf.dat.list) <- genomes
  return(gtf.dat.list)
}

################################################################################
#                          Prep Experimental Table                             #
################################################################################

# Experiment table will be iterated through
PrepExpTab <- function(input.file, star.dir, id.col){
  # Read in file
  xtab <- read.csv(input.file)
  # Add bam file locs to table
  xtab$bam <- file.path(star.dir, paste0(xtab[[id.col]], '.Aligned.sortedByCoord.out.bam'))
  # Add splice file locs to table
  xtab$splice <- file.path(star.dir, paste0(xtab[[id.col]], '.SJ.out.tab'))
  return(xtab)
}

################################################################################
#                                 Process BAM                                  #
################################################################################

ProcessBAMs_v1 <- function(obj.dir, xtab, id.col, gtf.list) {
  
  # Define output files
  xtab_rdat.file <- file.path(obj.dir, "virus_xtab_v1.rdat")
  junc_rdat.file <- file.path(obj.dir, "virus_GenomicAlignments_junction_discovery_v1.rdat")
  segment_rdat.file <- file.path(obj.dir, "virus_GenomicAlignments_segments_count_v1.rdat")
  
  # Output columns
  xtab$total_reads <- 0
  xtab$spliced_reads <- 0
  xtab$positive_strand <- 0
  xtab$negative_strand <- 0
  xtab$genome_reads <- 0
  
  # Only run if files don't already exist
  if (file.exists(xtab_rdat.file) &
      file.exists(junc_rdat.file) & 
      file.exists(segment_rdat.file)) {
    
    load(xtab_rdat.file)
    load(junc_rdat.file)
    load(segment_rdat.file)
    
  } else {
    
    # Prep empty lists
    segment_cnt.list <- list()
    junc.list <- list()
    
    # Loop through all samples
    for (r in 1:nrow(xtab)) {
      
      # Log
      message(paste('Starting iteration', r, 'of', nrow(xtab)))
      
      # Read in just the second mate and reverse the read-strand
      message('Reading BAM...')
      single.ga <- readGAlignments(
        file = BamFile(xtab$bam[r]),
        use.names = TRUE,
        param = ScanBamParam(
          what = c("seq"),
          #which = gtf.list[[xtab$virus[r]]],
          flag = scanBamFlag(
            isSecondMateRead = NA#,
            #isSecondaryAlignment = FALSE,
            #isSupplementaryAlignment = FALSE
          )
        ),
        with.which_label = FALSE
      )
      
      # Count total number of reads
      message('Counting total number of reads...')
      xtab$total_reads[r] <- length(single.ga)

      # Count number of +ve and -ve strand reads
      message('Counting +ve & -ve reads...')
      strand_info <- strand(single.ga)
      xtab$positive_strand[r] <- sum(strand_info == "+")
      xtab$negative_strand[r] <- sum(strand_info == "-")
      
      # Limit by strand
      single.ga <- single.ga[strand(single.ga) %in% "+", ]
      
      # Isolate spliced and unspliced
      message('Splitting spliced and unspliced...')
      spliced.ga <- single.ga[grepl("N", cigar(single.ga)), ]
      unspliced.ga <- single.ga[!grepl("N", cigar(single.ga)), ]
      
      # Count number of spliced and unspliced reads
      message('Counting spliced and unspliced...')
      xtab$spliced_reads[r] <- length(spliced.ga)
      xtab$unspliced_reads[r] <- length(unspliced.ga)
      
      # Count number of genome reads (for background)
      message('Counting genome reads...')
      # Define the genome background region as +/-500bp from the start of the S sgRNA intron
      S.Start <- max(gtf.list[[xtab$virus[r]]][gtf.list[[xtab$virus[r]]]$gene_id == 'S']@ranges@start)
      ORF1ab.Start <- gtf.list[[xtab$virus[r]]][gtf.list[[xtab$virus[r]]]$gene_id == 'ORF1ab' &
                                                  gtf.list[[xtab$virus[r]]]$type == 'exon']@ranges@start
      midpoint <- round((S.Start + ORF1ab.Start) / 2)
      genomeRegion <- GRanges(gtf.list[[xtab$virus[r]]]@seqnames@values,
                              IRanges(start = midpoint, width = 1000),
                              strand = "+")
      # Count unspliced reads that overlap
      genome_bg.ga <- single.ga[single.ga %over% genomeRegion &
                                  !grepl("N", cigar(single.ga)), ]
      xtab$genome_reads[r] <- length(genome_bg.ga)
      
      # Segment counts
      message('Calculating segment counts...')
      total.spliced.cnt <- table(seqnames(single.ga))[seqlevels(gtf.list[[xtab$virus[r]]])]
      seg.spliced.cnt <- table(seqnames(spliced.ga))[seqlevels(gtf.list[[xtab$virus[r]]])]
      seg.df <- data.frame(
        segment = seqlevels(gtf.list[[xtab$virus[r]]]),
        total   = as.numeric(total.spliced.cnt),
        spliced = as.numeric(seg.spliced.cnt)
      )
      segment_cnt.list[[as.character(xtab[[id.col]])[r]]] <- seg.df
      
      # Summarise junctions
      sample_junc.gr <- summarizeJunctions(single.ga, with.revmap = FALSE)
      junc.list[[as.character(xtab[[id.col]])[r]]] <- sample_junc.gr
      
      # Log
      message(paste('Finished iteration', r, 'of', nrow(xtab)))
      
    }
    
    save(xtab, file = xtab_rdat.file)
    save(junc.list, file = junc_rdat.file)
    save(segment_cnt.list, file = segment_rdat.file)
    
  }
}

ProcessBAMs_v2 <- function(obj.dir, xtab, id.col, gtf.list) {
  # Define output files
  xtab_rdat.file <- file.path(obj.dir, "virus_xtab_v2.rdat")
  junc_rdat.file <- file.path(obj.dir,
                              "virus_GenomicAlignments_junction_discovery_v2.rdat")
  segment_rdat.file <- file.path(obj.dir, "virus_GenomicAlignments_segments_count_v2.rdat")
  # Output columns
  xtab$total_reads <- 0
  xtab$spliced_reads <- 0
  xtab$positive_strand <- 0
  xtab$negative_strand <- 0
  xtab$genome_reads <- 0
  # Prep empty lists
  segment_cnt.list <- list()
  junc.list <- list()
  # Loop through all samples
  for (r in 1:nrow(xtab)) {
    # Log
    message(paste('Starting iteration', r, 'of', nrow(xtab)))
    # Create a sample_seg_cnt_list to add seg counts from each chunk too
    sample_seg_cnt_list <- list()
    # Create a sample_junc_list to add juncs from each chunk too
    sample_junc_list <- list()
    # Load bam in chunks
    chunkCnt <- 1
    bam <- BamFile(xtab$bam[r], yieldSize = 1e7)
    open(bam)
    repeat {
      message(paste('Starting chunk', chunkCnt))
      # Load chunk
      chunk <- readGAlignments(
        bam,
        use.names = FALSE,
        param = ScanBamParam(
          what = c("rname", "strand", "pos", "cigar", "flag", "qwidth"),
          which = gtf.list[[xtab$virus[r]]],
          flag = scanBamFlag(isSecondMateRead = NA)
        ),
        with.which_label = FALSE
      )
      # Break if chunk is empty
      if (length(chunk) == 0L)
        break
      # Process chunk here
      # Count number of +ve and -ve strand reads
      strand_info <- strand(chunk)
      xtab$positive_strand[r] <- xtab$positive_strand[r] + sum(strand_info == "+")
      xtab$negative_strand[r] <- xtab$negative_strand[r] + sum(strand_info == "-")
      # Limit by strand
      single.ga <- chunk[strand(chunk) %in% "+", ]
      # Isolate spliced and unspliced
      message('Splitting spliced and unspliced...')
      spliced.ga <- single.ga[grepl("N", cigar(single.ga)), ]
      unspliced.ga <- single.ga[!grepl("N", cigar(single.ga)), ]
      # Count number of reads in each
      message('Counting spliced and unspliced...')
      xtab$total_reads[r] <- xtab$total_reads[r] + length(single.ga)
      xtab$spliced_reads[r] <- xtab$spliced_reads[r] + length(spliced.ga)
      # Count number of genome reads (for background)
      message('Counting genome reads...')
      # Define the genome background region as +/-500bp from the start of the S sgRNA intron
      S.Start <- max(gtf.list[[xtab$virus[r]]][gtf.list[[xtab$virus[r]]]$gene_id == 'S']@ranges@start)
      ORF1ab.Start <- gtf.list[[xtab$virus[r]]][gtf.list[[xtab$virus[r]]]$gene_id == 'ORF1ab' &
                                                  gtf.list[[xtab$virus[r]]]$type == 'exon']@ranges@start
      midpoint <- round((S.Start + ORF1ab.Start) / 2)
      genomeRegion <- GRanges(gtf.list[[xtab$virus[r]]]@seqnames@values,
                              IRanges(start = midpoint, width = 1000),
                              strand = "+")
      # Count unspliced reads that overlap
      genome_bg.ga <- single.ga[single.ga %over% genomeRegion &
                                 !grepl("N", cigar(single.ga)), ]
      xtab$genome_reads[r] <-  xtab$genome_reads[r] + length(genome_bg.ga)
      # Segment counts
      message('Calculating segment counts...')
      total.spliced.cnt <- table(seqnames(single.ga))[seqlevels(gtf.list[[xtab$virus[r]]])]
      seg.spliced.cnt <- table(seqnames(spliced.ga))[seqlevels(gtf.list[[xtab$virus[r]]])]
      seg.df <- data.frame(
        segment = seqlevels(gtf.list[[xtab$virus[r]]]),
        total   = as.numeric(total.spliced.cnt),
        spliced = as.numeric(seg.spliced.cnt)
      )
      sample_seg_cnt_list[[length(sample_seg_cnt_list) + 1]] <- seg.df
      # Summarise junctions
      chunk_junc.gr <- summarizeJunctions(single.ga, with.revmap = FALSE)
      sample_junc_list[[length(sample_junc_list) + 1]] <- chunk_junc.gr
      # Iterate chunk count
      chunkCnt <- chunkCnt + 1
    }
    close(bam)
    # Summarise the segment counts for sample
    summary_df <- bind_rows(sample_seg_cnt_list) %>%
      group_by(segment) %>%
      summarise(
        total = sum(total, na.rm = TRUE),
        spliced = sum(spliced, na.rm = TRUE),
        .groups = "drop"
      )
    segment_cnt.list[[as.character(xtab[[id.col]])[r]]] <- summary_df
    # Summarise the junc list for sample
    junc_all <- do.call(c, sample_junc_list)
    # Ensure metadata column containing counts exists; common names: score, count, njunc
    # Replace "score" below with the correct metadata column name (e.g., "score" or "njunc")
    meta_name <- "score"
    stopifnot(meta_name %in% colnames(mcols(junc_all)))
    df <- as.data.frame(junc_all) %>%
      select(seqnames, start, end, strand, !!sym(meta_name)) %>%
      rename(value = !!sym(meta_name))
    summary_df <- df %>%
      group_by(seqnames, start, end, strand) %>%
      summarise(value = sum(value, na.rm = TRUE), .groups = "drop")
    # Rebuild a GRanges
    junc_merged <- GRanges(
      seqnames = summary_df$seqnames,
      ranges = IRanges(start = summary_df$start, end = summary_df$end),
      strand = summary_df$strand,
      score = summary_df$value
    )
    junc.list[[as.character(xtab[[id.col]])[r]]] <- junc_merged
    # Log
    message(paste('Finished iteration', r, 'of', nrow(xtab)))
  }
  save(xtab, file = xtab_rdat.file)
  save(junc.list, file = junc_rdat.file)
  save(segment_cnt.list, file = segment_rdat.file)
}


################################################################################
#                           Table of all junctions                             #
################################################################################

# Compile a table of all junctions with their depth of coverage.
JunctionTable <- function(xtab, id.col, junc.list) {
  
  # Empty list for downstream
  sample_junc.list <- list()
  
  # Loop through junction list
  for (r in 1:length(junc.list)) {
    # Log
    message(paste('Starting', r, 'of', length(junc.list)))
    
    sample_junc.gr <- junc.list[[r]]
    
    if (length(sample_junc.gr) == 0) {
      sample_junc.df <- data.frame(
        Sample_ID         = names(junc.list)[r],
        chr            = as.character(xtab$virus[match(names(junc.list)[r], as.character(xtab[[id.col]]))]),
        donor_site     = -1,
        acceptor_site  = -1,
        junction_depth = 0,
        stringsAsFactors = F
      )
      
    } else {
      sample_junc.df <- data.frame(
        Sample_ID         = names(junc.list)[r],
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
    
    sample_junc.df$total_junction_depth[is.na(sample_junc.df$total_junction_depth) |
                                          is.nan(sample_junc.df$total_junction_depth)] <- 0
    sample_junc.df$junction_freq[is.na(sample_junc.df$junction_freq) |
                                   is.nan(sample_junc.df$junction_freq)] <- 0
    sample_junc.df$log10_junction_freq[is.na(sample_junc.df$log10_junction_freq) |
                                         is.nan(sample_junc.df$log10_junction_freq)] <- 0
    
    sample_junc.list[[names(junc.list)[r]]] <- sample_junc.df
  }
  
  # Bind lists and dataframes
  junc.dat <- do.call(rbind, sample_junc.list)
  junc.dat <- junc.dat %>%
    left_join(xtab, by = join_by(Sample_ID))

  # Add gap size
  # Remove junction depth = 0
  junc.dat <- junc.dat %>%
    mutate(gap_size = acceptor_site - donor_site) %>%
    filter(junction_depth > 0)
  
  return(junc.dat)
  
}

################################################################################
#                      Add sequence of donor/acceptor                          #
################################################################################

AddJuncSeq <- function(junc.tab, fasta){
  
  # Index fasta
  indexFa(fasta)
  
  # Load FASTA file
  fastaObj <- FaFile(fasta)

  # Function to get donor/acceptor sequences for a given subset of junc.dat
  get_junc_seqs <- function(subdat, fa) {
    # Donor GRanges
    donor.gr <- GRanges(
      seqnames = subdat$chr,
      ranges   = IRanges(start = subdat$donor_site - 5, 
                         end   = subdat$donor_site + 20),
      strand = "+"
    )
    
    # Acceptor GRanges
    acceptor.gr <- GRanges(
      seqnames = subdat$chr,
      ranges   = IRanges(start = subdat$acceptor_site - 5, 
                         end   = subdat$acceptor_site + 20),
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
  
  # Save
  saveRDS(junc.seq.dat, file = paste0(obj.dir, '/virus_junction_seq_dat.RDS'))
  write.csv(junc.seq.dat, file = paste0(res.dir, '/virus_junction_seq_dat.csv'))
  
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

JuncClass <- function(xtab, id.col, gtf.list, junc.seq.dat) {
  
  # Define canonical junctions from the gtf file(s).
  canonical_junctions.gr <- list()
  for (n in 1:nrow(xtab)) {
    virus <- as.character(xtab$virus[n])
    gtf.gr <- gtf.list[[virus]]

    ## Canonical
    canonical.gr <- lapply(split(gtf.gr, gtf.gr$gene_id), function(x) {
      x[x$type %in% "exon", ]
    })
    canonical.gr <- canonical.gr[sapply(canonical.gr, length) > 1]      # Only splicing sgRNA
    for (g in 1:length(canonical.gr)) {
      canonical.tx.gr <- split(canonical.gr[[g]], canonical.gr[[g]]$transcript_id)
      canonical.tx.gr <- canonical.tx.gr[lapply(canonical.tx.gr, length) >
                                           1]
      for (t in 1:length(canonical.tx.gr)) {
        canonical.tx.gaps.gr <- gaps(canonical.tx.gr[[t]])
        canonical.tx.gaps.gr$virus <- virus
        canonical.tx.gaps.gr$gene_id <- names(canonical.gr)[g]
        canonical.tx.gaps.gr$transcript_name <- names(canonical.tx.gr)[t]
        canonical_junctions.gr <- c(canonical_junctions.gr, canonical.tx.gaps.gr)
      }
    }
  }
  canonical_junctions.gr <- unlist(GRangesList(canonical_junctions.gr))
  canonical_junctions.gr$key <- paste(
    canonical_junctions.gr$virus,
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
    canonical.idx <- canonical_junctions.gr$virus %in% junc.seq.dat$virus[r] &
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
      alternative.idx <- canonical_junctions.gr$virus %in% junc.seq.dat$virus[r] &
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
  
  # Save
  write.csv(junc.seq.dat,
            file = paste0(res.dir, '/virus_junction_seq_class.csv'))
  
  # Class counts weighted by read-depth
  class.tab <- t(as.data.frame.matrix(sapply(split(junc.seq.dat, junc.seq.dat[[id.col]]), function(x) {
    table(rep(x$class, x$junction_depth))
  })))
  perc.tab <- 100 * class.tab / rowSums(class.tab, na.rm = T)
  colnames(perc.tab) <- paste(colnames(perc.tab), ".perc", sep = '')
  perc.tab[is.nan(perc.tab)] <- 0
  class.df <- data.frame(
    sample_name  = rownames(class.tab),
    class.tab,
    total = rowSums(class.tab),
    perc.tab,
    stringsAsFactors = F
  )
  
  # Save
  write.csv(class.df,
            file = paste0(res.dir, '/virus_junction_class_summary.csv'))
  
  return(class.df)
  
}

################################################################################
#                       Junction Class Sequence Based                          #
################################################################################

## Define classes of spliced junction reads based on sequence motifs
# TRS-L = AATCTAAACTT
# NS2A = GTAATCTAAAC
# HE = AATATTAAACT
# S = TAATCTAAAC
# NS5A = GTAATCAAAACTT
# E = AATCCAAACATTAT
# M = GTAATCCAAACATTAT
# N = ATCTAAATTTTA

JuncClass_SeqBased <- function(junc.seq.dat) {
  
  TRSL = 'AATCTAAACTT'
  NS2A = 'GTAATCTAAAC'
  HE = 'AATATTAAACT'
  S = 'TAATCTAAAC'
  NS5A = 'GTAATCAAAACTT'
  E = 'AATCCAAACATTAT'
  M = 'GTAATCCAAACATTAT'
  N = 'ATCTAAATTTTA'
  
  sgRNAs <- c('NS2A', 'HE', 'S', 'NS5A', 'E', 'M', 'N')
  
  junc.class <- junc.seq.dat %>%
    mutate('Donor_Class' = case_when(grepl(TRSL, donor_seq) ~ 'TRS-L',
                                     .default = 'Unknown')) %>%
    mutate('Acceptor_Class' = case_when(grepl(NS2A, acceptor_seq) ~ 'NS2A',
                                        grepl(HE, acceptor_seq) ~ 'HE',
                                        grepl(S, acceptor_seq) ~ 'S',
                                        grepl(NS5A, acceptor_seq) ~ 'NS5A',
                                        grepl(E, acceptor_seq) ~ 'E',
                                        grepl(M, acceptor_seq) ~ 'M',
                                        grepl(N, acceptor_seq) ~ 'N',
                                        .default = 'Unknown')) %>%
    mutate('Class1' = case_when(Donor_Class == 'TRS-L' & Acceptor_Class %in% sgRNAs ~ 'Canonical',
                               Donor_Class == 'Unknown' & Acceptor_Class %in% sgRNAs ~ 'Alternative',
                               Donor_Class == 'TRS-L' & Acceptor_Class == 'Unknown' ~ 'Alternative',
                               Donor_Class == 'Unknown' & Acceptor_Class == 'Unknown' ~ 'Defective')) %>%
    mutate('Class2' = case_when(Class1 == 'Defective' & gap_size <= 10 ~ 'Defective_S1',
                                Class1 == 'Defective' & gap_size > 10 & gap_size <= 100 ~ 'Defective_S2',
                                Class1 == 'Defective' & gap_size > 100 & gap_size <= 1000 ~ 'Defective_S3',
                                Class1 == 'Defective' & gap_size > 1000 & gap_size <= 10000 ~ 'Defective_S4',
                                Class1 == 'Defective' & gap_size > 10000 & gap_size <= 100000 ~ 'Defective_S5',
                                Class1 == 'Canonical' ~ Acceptor_Class,
                                .default = Class1))
  
  # Save
  write.csv(junc.class,
            file = paste0(res.dir, '/virus_junction_seq_class.csv'))
  
  return(junc.class)
  
}
  

