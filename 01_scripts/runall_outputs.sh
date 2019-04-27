#!/bin/bash
# Run this script from within the stacks_workflow main repo

# 1. Identify loci out of HWE
echo "Identifying loci out of HWE"
# populations --in-path 05-stacks --popmap 01-info_files/population_map_retained.txt \
#         --out-path 07-filtered_vcfs -t 8 \
#         -r 0.7 --min-populations 15 --min-maf 0.01 \
#         --hwe 


# Make file of loci out of HWE
echo "Making list of markers out of HWE"
# awk -F"\t" ' $20 < 0.05 { print $1 "\t" $20 }' 07-filtered_vcfs/populations.sumstats.tsv  |\
#         grep -vE '^#' - > 07-filtered_vcfs/loci_out_of_hwe.txt
# awk -F"\t" '{ print $1 }' 07-filtered_vcfs/loci_out_of_hwe.txt | uniq > 07-filtered_vcfs/loci_out_of_hwe_blacklist.txt
echo "List name: 07-filtered_vcfs/loci_out_of_hwe_blacklist.txt ; use this as blacklist"

# 2. Filter and output loci for population genetic analysis (single SNP)
echo "Filtering and outputing single SNP per locus using blacklist"
# populations --in-path 05-stacks --popmap 01-info_files/population_map_retained.txt \
#         --out-path 07-filtered_vcfs -t 8 \
#         -r 0.7 --min-populations 15 --min-maf 0.01 \
#         --vcf --fstats --smooth --plink \
#         --write-single-snp \
#         --hwe \
#         -B 07-filtered_vcfs/loci_out_of_hwe_blacklist.txt

echo "Saving single-snp data"
# mkdir 09-diversity_stats 2>/dev/null 
# cp 07-filtered_vcfs/populations.snps.vcf 09-diversity_stats
# 
# mkdir 11-adegenet_analysis 2>/dev/null
# cp 07-filtered_vcfs/*plink* 11-adegenet_analysis 


# 3. Filter and output loci for population genetic analysis (single SNP)
echo "Filtering and outputting haplotypes per locus using blacklist"
populations --in-path 05-stacks --popmap 01-info_files/population_map_retained.txt \
        --out-path 07-filtered_vcfs -t 8 \
        -r 0.7 --min-populations 15 --min-maf 0.01 \
        --vcf --fstats --smooth --plink \
        --hwe \
        -B 07-filtered_vcfs/loci_out_of_hwe_blacklist.txt

echo "Saving haplotype and all data to 08-fineRADstructure"
mkdir 08-fineRADstructure 2>/dev/null
cp 07-filtered_vcfs/populations.haps.vcf 08-fineRADstructure

echo "Saving all data to 21_haplotype_results"
mkdir 21_haplotype_results 2>/dev/null
mv 07-filtered_vcfs/* 21_haplotype_results

# 4. Filter and output loci for selection analysis (single SNP)
echo "Filtering and outputting single SNP for selection analysis, no HWE filter"
populations --in-path 05-stacks --popmap 01-info_files/population_map_retained.txt \
        --out-path 07-filtered_vcfs -t 8 \
        -r 0.7 --min-populations 15 --min-maf 0.01 \
        --vcf --fstats --smooth --plink \
        --hwe 

echo "Saving all data to 24-selection"
mkdir 24-selection 2>/dev/null
mv 07-filtered_vcf/* 24-selection


echo "Done!"
