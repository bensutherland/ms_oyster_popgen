#!/bin/bash
# Run this script from within the stacks_workflow main repo
# Only works properly if there is only one pstacks.log file in the 10-log_files folder


# Obtain the filename
grep -E 'Alignments' 10-log_files/*_1b_pstacks.log |
    awk -F'/' '{print $2}' | sed 's/\.bam//g' - > temp_filenames.txt

# Obtain the number loci and mean coverage per sample
grep -E '^Kept' 10-log_files/*_1b_pstacks.log |
    awk -F' ' '{ print $2","$7 }' > temp_num_loci_and_mean_cov.txt

# Combine these two datasets
paste -d ',' temp_filenames.txt temp_num_loci_and_mean_cov.txt > temp_file.txt

# Add header
echo "sample,num.loci,avg.cov" > header.txt
cat header.txt temp_file.txt > output_pstacks_results.csv

# Clean workspace
rm temp_filenames.txt temp_num_loci_and_mean_cov.txt temp_file.txt header.txt

# Graph results
Rscript ./../ms_oyster_popgen/01_scripts/01b_assess_pstacks.R
