#!/bin/bash
# Set up and run bayescan

# Set variables
INFO_FILES_FOLDER="../ms_oyster_popgen/00_archive"
SET_CONTRASTS_FILE="selection_contrasts.txt"
SEL_FOLDER="13_selection"
SET_FOLDER="sets"
BAYESCAN_FORMAT_SCRIPT="/home/ben/Desktop/sfg_scripts/make_bayescan_input.py"
BAYESCAN="/home/ben/programs/BayeScan2.1/binaries/BayeScan2.1_linux64bits"

# Changing variables
THREADS=7


# Prepare inputs for selection
ls $SEL_FOLDER/*_strata.txt |
    while read f
    do

        # Identify the contrast
        CONTRAST=${f%_strata.txt}
        CONTRAST=$(basename $CONTRAST)
        echo "Working on the following contrast: $CONTRAST"

        # Create pop map files for each contrast
        echo "Creating pop_map from $f"
        CONTRAST_POP_MAP=${f%_strata.txt}"_pop_map.txt"
        cat $f | grep -vE '^INDIVIDUALS' > $CONTRAST_POP_MAP

        # Make bayescan input for each contrast
        echo "Make BayeScan input from VCF and pop map"
        CONTRAST_VCF=${f%_strata.txt}".recode.vcf"
        $BAYESCAN_FORMAT_SCRIPT $CONTRAST_VCF $CONTRAST_POP_MAP 10
        
        # Move the output
        mkdir $SEL_FOLDER/$CONTRAST
        mv bayes_input.txt $CONTRAST"_bayes_input.txt"
        mv $CONTRAST"_bayes_input.txt" low_freq_snps.txt snpkey.txt $SEL_FOLDER/$CONTRAST

        # Run BayeScan
        cd $SEL_FOLDER/$CONTRAST
        $BAYESCAN $CONTRAST"_bayes_input.txt" -o $CONTRAST -threads $THREADS -pr_odds 100
        cd ../..

        echo "Done working on $CONTRAST"
        
done
