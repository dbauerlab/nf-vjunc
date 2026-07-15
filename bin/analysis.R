#!/usr/bin/env Rscript

# Analysis script for quantification and classification of viral junctions/splicing/recombination 

# Script will work on a by sample basis unlike the original versions
# Will take the aligned sorted BAM file and the sample matches GTF as input and then output tables of junctions and segments for each sample which can then be combined for downstream analysis
# We will need to bind the sample BAM and GTF together in nextflow and then pass them together to this script for processing

# Libraries
library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(Biostrings)
library(Rsamtools)
library(reshape2)
library(stringr)
library(scales)
library(dplyr)
library(tidyr)

# Source custom functions
source("functions.R")

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
bam_file <- args[2]
gtf_file <- args[3]
fasta_file <- args[4]

# Import GTF
gtf <- rtracklayer::import(gtf_file)

# Process BAMs
results.list <- ProcessBAMs(sample = sample, bam = bam_file, gtf = gtf)

# Create junction table
juncTab <- JunctionTable(sample = sample, results.list = results.list)

# Add junction sequences
juncTabSeq <- AddJuncSeq(
  sample = sample,
  junc.tab = juncTab,
  fasta = fasta_file
)

# Add junction classes
juncTabClass <- JuncClass(sample = sample, junc.seq.dat = juncTabSeq, gtf = gtf)
