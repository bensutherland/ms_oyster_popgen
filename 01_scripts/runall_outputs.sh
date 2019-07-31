#!/bin/bash
# Run this script from within the stacks_workflow main repo
# A blacklist will be made using the populations HWE filter, moved to: 01-info_files/loci_out_of_hwe_blacklist.txt

# User set variables
POP_MAP="01-info_files/population_map_retained_names.txt"


#### 1. Identify loci out of HWE ####
#Identify markers out of HWE using the populations module 
echo "Running Populations to identify loci out of HWE"
populations --in-path 05-stacks --popmap $POP_MAP \
        --out-path 07-filtered_vcfs -t 8 \
        -r 0.7 --min-populations 15 --min-maf 0.01 \
        --hwe 

# Obtain header of sumstats file
echo "Collecting header for sumstats file"
grep -E '^#' 07-filtered_vcfs/populations.sumstats.tsv | tail -n 1 > 07-filtered_vcfs/populations.sumstats_header.txt

# Make blacklist with loci out of HWE
R --slave -e 'source("../ms_oyster_popgen/01_scripts/hwe_count_per_pop.R")'
echo "Saving the blacklist to 01-info_files"

# Identify most recent blacklist (to use)
NUM_BLACKLIST=$(ls -lth 01-info_files/blacklist* | awk -F"01-info_files/" '{ print $2 }' - | wc -l)
BLACKLIST=$(ls -lth 01-info_files/blacklist* | awk '{ print $9 }' - | head -n 1)

# Reporting on blacklist
if [ NUM_BLACKLIST > 1 ]
then
    echo "Warning, there are multiple blacklists identified, using the most recent"
    echo $BLACKLIST
else
    echo "Using the blacklist $BLACKLIST"
fi

#### 2. Filter and output loci for population genetic analysis (single SNP) ####
echo "Running Populations to output single SNP data (with blacklist)"
populations --in-path 05-stacks --popmap $POP_MAP \
        --out-path 07-filtered_vcfs -t 8 \
        -r 0.7 --min-populations 15 --min-maf 0.01 \
        --vcf --fstats --smooth --plink \
        --fasta-loci \
        --write-single-snp \
        --hwe \
        -B $BLACKLIST

echo "Saving single-snp PLINK data - Differentiation"
mkdir 11-adegenet_analysis 2>/dev/null
cp 07-filtered_vcfs/*plink* 11-adegenet_analysis 

echo "Saving all other single-snp data - Diversity"
mkdir 09-diversity_stats 2>/dev/null 
cp 07-filtered_vcfs/* 09-diversity_stats


#### 3. Filter and output loci for population genetic analysis (haplotypes) ####
echo "Running Populations to output haplotypes data (using blacklist)"
populations --in-path 05-stacks --popmap $POP_MAP \
        --out-path 07-filtered_vcfs -t 8 \
        -r 0.7 --min-populations 15 --min-maf 0.01 \
        --vcf --fstats --smooth --plink \
        --hwe \
        --radpainter \
        -B $BLACKLIST 

echo "Saving haplotypes to 08-fineRADstructure"
mkdir 08-fineRADstructure 2>/dev/null
cp 07-filtered_vcfs/populations.haps.vcf 08-fineRADstructure

echo "Saving all data to 21_haplotype_results"
mkdir 21-haplotype_results 2>/dev/null
mv 07-filtered_vcfs/* 21-haplotype_results
 
echo "Done all outputs!"
