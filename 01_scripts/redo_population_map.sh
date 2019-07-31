#!/bin/bash
# Identify what samples are truly in your study, and rebuild the population map using only these files. Following this, populations can be merged using rename_populations_in_map.R

# Define variables
INFO_DIR="01-info_files"
SAMPLE_DIR="04-all_samples"

# First, create a file that lists only the retained fastq files in 04-all_samples
ls $SAMPLE_DIR/*.fq.gz | awk -F"/" '{ print $2 }' - | sed 's/\.fq\.gz//g' - > $INFO_DIR/all_retained_samples.txt

# Second, keep only the lines of the population map that are in the retained samples
grep -f $INFO_DIR/all_retained_samples.txt $INFO_DIR/population_map.txt > $INFO_DIR/population_map_retained_numeric.txt

