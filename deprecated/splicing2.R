# Script for splicing analysis
# Adapted from Richard Mitters

################################################################################
#                                    Setup                                     #
################################################################################

# Packages
library(ggplot2)
library(kableExtra)
library(tximport)
library(limma)
library(DESeq2)
library(RColorBrewer)
library(viridis)
library(GenomicFeatures)
library(GenomicAlignments)
library(reshape2)
library(patchwork)
library(circlize)
library(magick)
library(stringr)
library(openxlsx)
library(scales)
library(DT)
library(Rsamtools)
library(Biostrings)
library(GenomicRanges)
library(dplyr)
library(tidyr)

# Fixed directories
proj.dir <- "/nemo/stp/babs/working/bootj/projects/bauerd/danyn.patel/flu_rnaseq"
data.dir <- paste(proj.dir, "/secondpass", sep = '')

# Output directories
obj.dir <- paste(proj.dir, "/objects", sep = '')
res.dir <- paste(proj.dir, "/results", sep = '')

# Create
if (!file.exists(obj.dir)) {
  dir.create(obj.dir)
}
if (!file.exists(res.dir)) {
  dir.create(res.dir)
}

# Find gtf files - virus only!
gtfFiles <- list.files(path = '../genomes',
                       pattern = "^B_.*\\.gtf$",
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

# Define the host genome name
host <- 'ROS_Cfam_1.0_'

################################################################################
#                          Prep Experimental Table                             #
################################################################################

# Experiment table will be iterated through
xtab <- read.csv('../samplesheets/experiment_table.csv')

# Add BAM files
xtab$bam <- paste0(proj.dir, "/secondpass/", xtab$species, '/results/star_salmon/', xtab$ID, '.markdup.sorted.bam')

# Add splice junctions
xtab$splice_juncs <- paste0(proj.dir, "/secondpass/", xtab$species, '/results/star_salmon/log/', xtab$ID, '.SJ.out.tab')

# Add quant.sf
xtab$salmon <- paste0(proj.dir, "/secondpass/", xtab$species, '/results/star_salmon/', xtab$ID, '/quant.sf')

# Add tx2gene
xtab$tx2gene <- paste0(proj.dir, "/secondpass/", xtab$species, '/results/star_salmon/tx2gene.tsv')

# Add virus
xtab$virus <- gsub('.gtf', '', xtab$gtf)

################################################################################
#                                 Process BAM                                  #
################################################################################

# Define output files 
xtab_rdat.file <- paste(obj.dir, "virus_xtab.rdat", sep = '/')
junc_rdat.file <- paste(obj.dir,
                            "virus_GenomicAlignments_junction_discovery.rdat",
                            sep = '/')
segment_rdat.file <- paste(obj.dir, "virus_GenomicAlignments_segments_count.rdat", sep =
                             '/')

# Output columns
xtab$total_reads <- 0
xtab$spliced_reads <- 0

# Only run if files don't already exist 
if (file.exists(xtab_rdat.file) & file.exists(junc_rdat.file) & file.exists(segment_rdat.file)) {
  
  load(xtab_rdat.file)
  load(junc_rdat.file)
  load(segment_rdat.file)
  
} else {
  
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
      index = file,
      use.names = FALSE,
      param = ScanBamParam(
        what = c("seq"),
        which = gtf.dat.list[[xtab$virus[r]]],
        flag = scanBamFlag(isSecondMateRead = TRUE)
      ),
      with.which_label = FALSE
    )
    # Limit by strand
    single.ga <- single.ga[strand(single.ga) %in% "+", ]
    
    # Isolate spliced and unspliced
    message('Splitting spliced and unspliced...')
    spliced.ga <- single.ga[grepl("N", cigar(single.ga)), ]
    unspliced.ga <- single.ga[!grepl("N", cigar(single.ga)), ]
    
    # Count number of reads in each
    message('Counting spliced and unspliced...')
    xtab$total_reads[r] <- length(single.ga)
    xtab$spliced_reads[r] <- length(spliced.ga)
    
    # Nucleotide coverage
    message('Calculating coverage...')
    total.cov   <- coverage(single.ga)
    spliced.cov <- coverage(spliced.ga)
    
    # Segment counts
    message('Calculating segment counts...')
    total.spliced.cnt <- table(seqnames(single.ga))[seqlevels(gtf.dat.list[[xtab$virus[r]]])]
    seg.spliced.cnt <- table(seqnames(spliced.ga))[seqlevels(gtf.dat.list[[xtab$virus[r]]])]
    seg.df <- data.frame(
      segment = seqlevels(gtf.dat.list[[xtab$virus[r]]]),
      total   = as.numeric(total.spliced.cnt),
      spliced = as.numeric(seg.spliced.cnt)
    )
    segment_cnt.list[[as.character(xtab$Sample)[r]]] <- seg.df
    
    # Summarise junctions
    sample_junc.gr <- summarizeJunctions(single.ga, with.revmap = FALSE)
    junc.list[[as.character(xtab$Sample)[r]]] <- sample_junc.gr
    
    # Log
    message(paste('Finished iteration', r, 'of', nrow(xtab)))
    
  }
  
  save(xtab,file=xtab_rdat.file)
  save(junc.list,file=junc_rdat.file)
  save(segment_cnt.list,file=segment_rdat.file)
  
}

################################################################################
#                           Table of all junctions                             #
################################################################################

# Compile a table of all junctions with their depth of coverage.

sample_junc.list <- list()

# Loop through junction list 
for (r in 1:length(junc.list)) {
  
  # Log
  message(paste('Starting', r, 'of', length(junc.list)))
  
  sample_junc.gr <- junc.list[[r]]
  
  if (length(sample_junc.gr) == 0) {
    sample_junc.df <- data.frame(
      sample         = names(junc.list)[r],
      chr            = as.character(xtab$virus[match(names(junc.list)[r], as.character(xtab$Sample))]),       
      donor_site     = -1,
      acceptor_site  = -1,
      junction_depth = 0,
      stringsAsFactors=F)
    
  } else {
    
    sample_junc.df <- data.frame(
      sample         = names(junc.list)[r],
      chr            = seqnames(sample_junc.gr),
      donor_site     = start(sample_junc.gr),
      acceptor_site  = end(sample_junc.gr),
      junction_depth = sample_junc.gr$score,
      stringsAsFactors=F)
  
    }
  
  # Junction frequency = junction depth / total junction depth
  sample_junc.df$total_junction_depth = sum(sample_junc.df$junction_depth)
  sample_junc.df$junction_freq        = sample_junc.df$junction_depth / sample_junc.df$total_junction_depth
  sample_junc.df$log10_junction_freq  = log10(sample_junc.df$junction_freq)
  
  sample_junc.df$total_junction_depth[is.na(sample_junc.df$total_junction_depth) | is.nan(sample_junc.df$total_junction_depth)] <- 0
  sample_junc.df$junction_freq[is.na(sample_junc.df$junction_freq) | is.nan(sample_junc.df$junction_freq)] <- 0
  sample_junc.df$log10_junction_freq[is.na(sample_junc.df$log10_junction_freq) | is.nan(sample_junc.df$log10_junction_freq)] <- 0
  
  sample_junc.list[[names(junc.list)[r]]] <- sample_junc.df
}

# Bind lists and dataframes
junc.dat <- do.call(rbind, sample_junc.list)
junc.dat <- cbind(junc.dat, xtab[match(junc.dat$sample, as.character(xtab$Sample)), c('total_reads', 'spliced_reads', 'virus')])

# Add gap size
# Remove junction depth = 0
junc.dat <- junc.dat %>%
  mutate(gap_size = acceptor_site - donor_site) %>%
  filter(junction_depth > 0)

################################################################################
#                      Add sequence of donor/acceptor                          #
################################################################################

# 1. Load all FASTA files into a named list of FaFile objects
fasta_files <- list.files(
  path = '../genomes',
  pattern = "^B_.*\\.fasta$",
  recursive = TRUE,
  full.names = TRUE
)

fastaList <- lapply(fasta_files, FaFile)
names(fastaList) <- gsub("\\.fasta$", "", basename(fasta_files))

# 2. Function to get donor/acceptor sequences for a given subset of junc.dat
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

junc.seq.dat <- do.call(
  rbind,
  lapply(split(junc.dat, junc.dat$virus), function(df) {
    fa <- fastaList[[df$virus[1]]]
    get_junc_seqs(df, fa)
  })
)

# Sort by junction freq
junc.seq.dat <- junc.seq.dat %>%
  arrange(desc(log10_junction_freq))

# Save
saveRDS(junc.seq.dat, file = paste0(obj.dir, '/virus_junction_seq_dat.RDS'))
write.csv(junc.seq.dat, file = paste0(res.dir, '/virus_junction_seq_dat.csv'))

################################################################################
#                               Gap Frequency                                  #
################################################################################

# Create a group variable 
junc.seq.dat$group <- gsub('_br[0-9]', '', junc.seq.dat$sample)

# Plot size distribution
plt <- ggplot(junc.seq.dat, aes(x = gap_size)) +
  geom_histogram(binwidth = 100, color = 'black', fill = 'lightblue') +
  theme_bw() +
  facet_grid(cols = vars(chr), rows = vars(group)) +
  theme(strip.text.x = element_text(size = 6, colour = "black", angle = 45),
        strip.text.y = element_text(size = 6, colour = "black", angle = 45),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
ggsave(plot = plt, filename = paste0(res.dir, '/gapsizes_hist.png'), height = 15, width = 15, units = 'in')

################################################################################
#                             Junction Class                                   #
################################################################################

## Define classes of spliced junction reads based on overlap with canonical mRNA.
## The junction reads are redefined as the "gaps" (introns) between the exons, to match the gtf definitions.
## canonical = matches 2/2 junctions
## alternative = matches 1/2 junctions
## defective = matches 0/2 junctions

## Define canonical junctions from the gtf file(s).
canonical_junctions.gr <- list()
for (n in 1:nrow(xtab)) {
  virus <- as.character(xtab$virus[n])
  gtf.gr <- gtf.dat.list[[virus]]
  junc.gr <- junc.list[[as.character(xtab$Sample[n])]]
  
  ## Canonical
  canonical.gr <- lapply(split(gtf.gr,gtf.gr$gene_id),function(x){x[x$type %in% "exon",]})
  canonical.gr <- canonical.gr[sapply(canonical.gr,length)>1]      # Only splicing sgRNA
  for (g in 1:length(canonical.gr)) {
    canonical.tx.gr <- split(canonical.gr[[g]],canonical.gr[[g]]$transcript_id)
    canonical.tx.gr <- canonical.tx.gr[lapply(canonical.tx.gr,length)>1]
    for (t in 1:length(canonical.tx.gr)) {
      canonical.tx.gaps.gr <- gaps(canonical.tx.gr[[t]])
      canonical.tx.gaps.gr$virus <- virus
      canonical.tx.gaps.gr$gene_id <- names(canonical.gr)[g]
      canonical.tx.gaps.gr$transcript_name <- names(canonical.tx.gr)[t]
      canonical_junctions.gr <- c(canonical_junctions.gr,canonical.tx.gaps.gr)
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
canonical_junctions.gr <- canonical_junctions.gr[!duplicated(canonical_junctions.gr$key),]

## Classify junctions

# Default is defective
junc.seq.dat$class <- factor("defective", levels = c("canonical", "alternative", "defective"))
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
      (((start(canonical_junctions.gr) %in% junc.seq.dat$donor_site[r]) &
          (!end(canonical_junctions.gr) %in% junc.seq.dat$acceptor_site[r])
      ) |
        ((!start(canonical_junctions.gr) %in% junc.seq.dat$donor_site[r]) &
           (end(canonical_junctions.gr) %in% junc.seq.dat$acceptor_site[r])
        ))
    
    if (any(alternative.idx)) {
      junc.seq.dat$class[r] <- "alternative"
      junc.seq.dat$gene_id[r] <- paste(unique(canonical_junctions.gr[alternative.idx,]$gene_id),collapse=';')
      junc.seq.dat$transcript_name[r] <- paste(unique(canonical_junctions.gr[alternative.idx,]$transcript_name),collapse=';')
    }
  }
}

# Save
write.csv(junc.seq.dat, file = paste0(res.dir, '/virus_junction_seq_class.csv'))

# Class counts weighted by read-depth
class.tab <- t(as.data.frame.matrix(sapply(split(junc.seq.dat, junc.seq.dat$sample), function(x) {
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
write.csv(class.df, file = paste0(res.dir, '/virus_junction_class_summary.csv'))

################################################################################
#                                 Plotting                                     #
################################################################################

# Plot n spliced reads
plt_tab <- junc.seq.dat %>%
  mutate(unspliced_reads = total_reads - spliced_reads) %>%
  select(sample, group, spliced_reads, unspliced_reads) %>%
  pivot_longer(!c(sample, group), names_to = 'reads', values_to = 'counts')

# Grouped bar
plt <- ggplot(plt_tab, aes(fill = reads, y = counts, x = sample)) +
  geom_bar(position = "dodge", stat = "identity") +
  facet_wrap(~ group, scales = "free_x", nrow = 1) +
  scale_y_log10(label = comma) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 4, angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size = 4),
        axis.title.x = element_text(size = 6),
        axis.title.y = element_text(size = 6),
        strip.text.x = element_text(size = 4, colour = "black", angle = 45),
        legend.text = element_text(size = 4), 
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.2, "cm")) +
  xlab('Sample') +
  ylab('log10(n reads)') + 
  guides(fill = guide_legend(title = "Read type"))
ggsave(plot = plt, filename = paste0(res.dir, '/nspliced_bar.png'), height = 3, width = 9, units = 'in')

# Plot splice classes across samples
class.df$group <- gsub('_br[0-9]', '', class.df$sample_name)

# Pivot
plt_tab <- class.df %>%
  select(sample_name, group, canonical, alternative, defective) %>%
  pivot_longer(!c(sample_name, group), names_to = 'Class', values_to = 'n_reads')

# Grouped bar
plt <- ggplot(plt_tab, aes(fill = Class, y = n_reads, x = sample_name)) +
  geom_bar(position = "dodge", stat = "identity") +
  facet_wrap(~ group, scales = "free_x", nrow = 1) +
  scale_y_log10(label = comma) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 4, angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size = 4),
        axis.title.x = element_text(size = 6),
        axis.title.y = element_text(size = 6),
        strip.text.x = element_text(size = 4, colour = "black", angle = 45),
        legend.text = element_text(size = 4), 
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.2, "cm")) +
  xlab('Sample') +
  ylab('log10(n_reads)') + 
  guides(fill = guide_legend(title = "Junction Class"))
ggsave(plot = plt, filename = paste0(res.dir, '/junctionClass_bar.png'), height = 3, width = 9, units = 'in')



