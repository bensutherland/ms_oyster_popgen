#!/bin/bash

# Determine number of reads per file from samples file
for i in $(ls 04-all_samples/*.fq.gz) ; do echo $i ; gunzip -c $i | grep -cE '^@' - ; done > reads_per_sample.txt

# Separate by second line into two columns
sed 'N;s/\n/ /' reads_per_sample.txt > reads_per_sample_table.txt

cat reads_per_sample_table.txt

# Run Rscript to calculate
#awk '{ print $2 }' reads_per_sample_table.txt | R --vanilla --slave -e "data=read.delim(pipe('cat /dev/stdin'), header =F) ; mean(data\$V1)"

# Better, just use this:
Rscript ./summary_stats.R

# How many alignments?



