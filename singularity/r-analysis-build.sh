#!/bin/bash

# Bash script for building the R environment for the analysis of viral junctions/splicing/recombination in nf-vjunc
# This script will be used in the singularity container to set up the R environment with the necessary packages for the analysis
# Note access token for singularity hub is required to build the container
# Mine is found: /camp/home/bootj/.singularity/remote.yaml

# Load singularity module
ml Singularity/3.11.3

# Disable bind paths for remote build
export SINGULARITY_BIND=""

# Build remotely
singularity build --remote r-analysis.sif r-analysis.def