#!/bin/bash
# Run this script from within the stacks_workflow main repo

# Set the samples to use for the comparisons
# Here using Hisnit, Pendrell, Pipestem and Serpentine for the BC wild
vcf-query -l 09-diversity_stats/batch_1.vcf | grep -E 'Hisnit|Pendrell|Pipestem|Serpentine' - | grep -vE 'PendrellFarm' - > 09-diversity_stats/wild_bc_samples.txt

# Here using the BC wild-to-farm sample
vcf-query -l 09-diversity_stats/batch_1.vcf | grep -E 'PendrellFarm' - > 09-diversity_stats/bc_wild_to_farm_samples.txt

# Here using Rosewall Farm
vcf-query -l 09-diversity_stats/batch_1.vcf | grep -E 'Rosewall' - > 09-diversity_stats/rosewall_samples.txt 

# Here for France
vcf-query -l 09-diversity_stats/batch_1.vcf | grep -E 'FranceW' - > 09-diversity_stats/france_wild_samples.txt
vcf-query -l 09-diversity_stats/batch_1.vcf | grep -E 'FranceC' - > 09-diversity_stats/france_cultured_samples.txt

# Here for China (altogether)
vcf-query -l 09-diversity_stats/batch_1.vcf | grep -E 'China' - > 09-diversity_stats/china_samples.txt

# Here for Guernsey
vcf-query -l 09-diversity_stats/batch_1.vcf | grep -E 'Guernsey' - > 09-diversity_stats/guernsey_samples.txt


# And Loop
ls -1 09-diversity_stats/*_samples.txt |
        perl -pe 's/\.txt//' |
        while read i
        do
          echo $i
          name=$(basename $i)
          vcftools --vcf 09-diversity_stats/batch_1.vcf --site-pi --keep $i".txt"
          mv out.sites.pi 09-diversity_stats/$name"_out.sites.pi"
done
