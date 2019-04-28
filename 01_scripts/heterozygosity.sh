#!/bin/bash
# Calculate heterozygosity with vcftools

DIVERSITY_FOLDER="09-diversity_stats"
VCF="populations.snps.vcf"


vcftools --vcf $DIVERSITY_FOLDER/$VCF --het    
# vcftools outputs the data to the working directory
mv out.het 09-diversity_stats/populations.snps.het
