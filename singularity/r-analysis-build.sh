#!/bin/bash

# Bash script for building the R environment for the analysis of viral junctions/splicing/recombination in nf-vjunc
# This script will be used in the singularity container to set up the R environment with the necessary packages for the analysis

ml Singularity/3.11.3

singularity build --fakeroot r-analysis.sif r-analysis.def