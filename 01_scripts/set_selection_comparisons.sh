#!/bin/bash
# Make restricted VCF files for contrasts of interest as specified using the SET_CONTRASTS_FILE

# Set variables
INFO_FILES_FOLDER="../ms_oyster_popgen/00_archive"
SET_CONTRASTS_FILE="selection_contrasts.txt"
SEL_FOLDER="13_selection"
SET_FOLDER="sets"
IN_VCF="populations.snps.vcf"


# Prepare inputs for selection
cat $INFO_FILES_FOLDER/$SET_CONTRASTS_FILE |
    while read f
    do
        # Define the two pops per contrast
        POP1=$(echo $f | awk '{ print $1 }' - ) 
        POP2=$(echo $f | awk '{ print $2 }' - ) 
        echo "Setting up contrast $POP1 vs $POP2" 

        SET_FILENAME="$POP1"_"$POP2".txt         

        # Create set lists
        grep -E "$POP1|$POP2" $SEL_FOLDER/$SET_FOLDER/all_samples.txt \
               > $SEL_FOLDER/$SET_FOLDER/$SET_FILENAME

        # Make restricted VCF files for each contrast
        vcftools --vcf $SEL_FOLDER/$IN_VCF \
            --keep $SEL_FOLDER/$SET_FOLDER/$SET_FILENAME \
            --out $SEL_FOLDER/${SET_FILENAME%.txt} \
            --recode
        # Cleanup
        mv $SEL_FOLDER/*.log $SEL_FOLDER/$SET_FOLDER/

done




