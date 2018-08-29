#!/bin/bash
# First obtain a list of samples containing keywords (part 1) then calculate Pi for each of these groups of samples.    
# Run this script from within the stacks_workflow main repo

DIVERSITY_FOLDER="09-diversity_stats"

# Part 1: Collect sample lists for comparison groups 
# BC Wild = Hisnit, Pendrell, Pipestem and Serpentine
vcf-query -l $DIVERSITY_FOLDER/batch_1.vcf |
    grep -E 'HIS_|PEN_|PIP_|SER_' - \
        > $DIVERSITY_FOLDER/bc_wild_samples.txt

# BC Wild-to-farm = PendrellFarm 
vcf-query -l $DIVERSITY_FOLDER/batch_1.vcf | 
    grep -E 'PENF_' - \
        > $DIVERSITY_FOLDER/bc_wild_to_farm_samples.txt

# BC Farm 1 = Rosewall 
vcf-query -l $DIVERSITY_FOLDER/batch_1.vcf | 
    grep -E 'ROS_' - \
        > $DIVERSITY_FOLDER/rosewall_samples.txt 

# BC Farm 2 = DeepBay
vcf-query -l $DIVERSITY_FOLDER/batch_1.vcf | 
    grep -E 'DPB_' - \
        > $DIVERSITY_FOLDER/deepbay_samples.txt 

# France Wild = FranceW
vcf-query -l $DIVERSITY_FOLDER/batch_1.vcf | 
    grep -E 'FRA_' - \
        > $DIVERSITY_FOLDER/france_wild_samples.txt

# France Wild-to-farm
vcf-query -l $DIVERSITY_FOLDER/batch_1.vcf | 
    grep -E 'FRAF_' - \
        > $DIVERSITY_FOLDER/france_wild_to_farm_samples.txt

# China Farm 1
vcf-query -l $DIVERSITY_FOLDER/batch_1.vcf | 
    grep -E 'QDC_' - \
        > $DIVERSITY_FOLDER/chinaQDC_samples.txt

# China Farm 2
vcf-query -l $DIVERSITY_FOLDER/batch_1.vcf | 
    grep -E 'RSC_' - \
        > $DIVERSITY_FOLDER/chinaRSC_samples.txt

# China Wild 
vcf-query -l $DIVERSITY_FOLDER/batch_1.vcf | 
    grep -E 'CHN_' - \
        > $DIVERSITY_FOLDER/china_wild_samples.txt

# China wild-to-farm
vcf-query -l $DIVERSITY_FOLDER/batch_1.vcf | 
    grep -E 'CHNF_' - \
        > $DIVERSITY_FOLDER/china_wild_to_farm_samples.txt

# Guernsey Farm
vcf-query -l $DIVERSITY_FOLDER/batch_1.vcf | 
    grep -E 'GUR_' - \
        > $DIVERSITY_FOLDER/guernsey_samples.txt

# Japan
vcf-query -l $DIVERSITY_FOLDER/batch_1.vcf | 
    grep -E 'JPN_' - \
        > $DIVERSITY_FOLDER/japan_samples.txt


# Part 2. Loop across the different 'samples' list files and calculate site pi
ls -1 09-diversity_stats/*_samples.txt |
        perl -pe 's/\.txt//' |
        while read i
        do
          echo $i
          name=$(basename $i)
          vcftools --vcf 09-diversity_stats/batch_1.vcf --site-pi --keep $i".txt"
          mv out.sites.pi 09-diversity_stats/$name"_out.sites.pi"
done
