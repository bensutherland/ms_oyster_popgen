#!/bin/bash
# Run this script from within the stacks_workflow main repo
# A blacklist will be made using the populations HWE filter, moved to: 01-info_files/loci_out_of_hwe_blacklist.txt

#### 1. Identify loci out of HWE ####
# Identify markers out of HWE using the populations module 
echo "Running Populations to identify loci out of HWE"
populations --in-path 05-stacks --popmap 01-info_files/population_map_retained.txt \
        --out-path 07-filtered_vcfs -t 8 \
        -r 0.7 --min-populations 15 --min-maf 0.01 \
        --hwe 

# Make file of loci out of HWE
echo "Making list of markers out of HWE"
awk -F"\t" ' $20 < 0.05 { print $1 "\t" $20 }' 07-filtered_vcfs/populations.sumstats.tsv  |\
        grep -vE '^#' - > 07-filtered_vcfs/loci_out_of_hwe.txt
awk -F"\t" '{ print $1 }' 07-filtered_vcfs/loci_out_of_hwe.txt | uniq > 07-filtered_vcfs/loci_out_of_hwe_blacklist.txt

# Save blacklist to info files folder
echo "Saving the blacklist to 01-info_files"
cp 07-filtered_vcfs/loci_out_of_hwe_blacklist.txt 01-info_files/loci_out_of_hwe_blacklist.txt

echo "Blacklist: 01-info_files/loci_out_of_hwe_blacklist.txt"


#### 2. Filter and output loci for population genetic analysis (single SNP) ####
echo "Running Populations to output single SNP data (with blacklist)"
populations --in-path 05-stacks --popmap 01-info_files/population_map_retained.txt \
        --out-path 07-filtered_vcfs -t 8 \
        -r 0.7 --min-populations 15 --min-maf 0.01 \
        --vcf --fstats --smooth --plink \
        --fasta-loci \
        --write-single-snp \
        --hwe \
        -B 01-info_files/loci_out_of_hwe_blacklist.txt

echo "Saving single-snp PLINK data - Differentiation"
mkdir 11-adegenet_analysis 2>/dev/null
cp 07-filtered_vcfs/*plink* 11-adegenet_analysis 

echo "Saving all other single-snp data - Diversity"
mkdir 09-diversity_stats 2>/dev/null 
cp 07-filtered_vcfs/* 09-diversity_stats


#### 3. Filter and output loci for population genetic analysis (haplotypes) ####
echo "Running Populations to output haplotypes data (using blacklist)"
populations --in-path 05-stacks --popmap 01-info_files/population_map_retained.txt \
        --out-path 07-filtered_vcfs -t 8 \
        -r 0.7 --min-populations 15 --min-maf 0.01 \
        --vcf --fstats --smooth --plink \
        --hwe \
        --radpainter \
        -B 07-filtered_vcfs/loci_out_of_hwe_blacklist.txt

echo "Saving haplotypes to 08-fineRADstructure"
mkdir 08-fineRADstructure 2>/dev/null
cp 07-filtered_vcfs/populations.haps.vcf 08-fineRADstructure

echo "Saving all data to 21_haplotype_results"
mkdir 21-haplotype_results 2>/dev/null
mv 07-filtered_vcfs/* 21-haplotype_results
 
 echo "Done all outputs!"
