#!/bin/bash
# Run this script from within the stacks_workflow main repo

# Move into the directory containing files
cd 04-all_samples 

# Clear previous runs
rm reads_per_sample.txt 2> /dev/null
rm reads_per_sample_table.txt 2> /dev/null
rm mappings_per_sample.txt 2> /dev/null 
rm mappings_per_sample_table.txt 2> /dev/null

## Number Reads
# Determine number of reads per file from samples file
for i in $(ls *.fq.gz) ; 
    do echo $i ;
    gunzip -c $i | grep -cE '^@' - ;
    done >> reads_per_sample.txt

# Separate by second line into two columns
sed 'N;s/\n/ /' reads_per_sample.txt > reads_per_sample_table.txt

# Remove intermediate file
rm reads_per_sample.txt


## Number Mappings
# Determine number of reads per file from samples file
for i in $(ls *.bam) ; 
    do echo $i ;
    samtools view $i | wc -l ;
    done >> mappings_per_sample.txt

# Separate by second line into two columns
sed 'N;s/\n/ /' mappings_per_sample.txt > mappings_per_sample_table.txt

# Remove intermediate file
rm mappings_per_sample.txt

# Run Rscript to calculate
#awk '{ print $2 }' reads_per_sample_table.txt | R --vanilla --slave -e "data=read.delim(pipe('cat /dev/stdin'), header =F) ; mean(data\$V1)"

cd ..

# Run some quick stats on the result 
Rscript ./../ms_oyster_popgen/01_scripts/summary_stats.R

