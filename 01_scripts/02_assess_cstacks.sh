#!/bin/bash
# This assesses the results of cstacks via the log file
# Run from main directory 

# Set variables
MY_FILE="cstacks_output_trimmed.txt"
TEMP_FILE="cstacks_output_trimmed_temp.txt"

# Obtain the necessary information and save as a temp file
grep -E 'newly added| were matched to a catalog locus\.' 10-log_files/*cstacks.log > $MY_FILE 

# Clean up input file
awk -F, '{if (NR!=1) { print $1} }' $MY_FILE > $TEMP_FILE 
# removes the unneeded piece after the comma, and the first sample, as this one doesn't match the rest of the file 

# Create temporary files for each column
awk 'NR % 2 == 1' $TEMP_FILE > col1.txt # number newly added 
awk 'NR % 2 == 0' $TEMP_FILE > col2.txt # number loci matched 
 
# Combine
paste -d ',' col1.txt col2.txt > $MY_FILE.table.txt 
 
rm col1.txt col2.txt $TEMP_FILE $MY_FILE

# Create output table
sed 's/ loci were matched to a catalog locus\.//g' $MY_FILE.table.txt |
    sed 's/ loci were newly added to the catalog\.//g' > cstacks_output_table.txt

echo "matched.to.catalog,newly.added" > header.csv

cat header.csv cstacks_output_table.txt > cstacks_output_table.csv

rm $MY_FILE.table.txt cstacks_output_table.txt header.csv
