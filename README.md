#### Analyze MOP oyster populations using Eric Normandeau's stacks_workflow pipeline
##### https://github.com/enormandeau/stacks_workflow

### Setup
Put all raw data into `02-raw` using cp or cp -l    
Make a new directory in the raw data folder
`mkdir 02-raw/fastqc_raw`    

Run fastqc and summarize results    
`fastqc 02-raw/*.fastq.gz -o 02-raw/fastqc_raw/ -t 5`    
`multiqc -o 02-raw/fastqc_raw/ 02-raw/fastqc_raw`   

Prepare lane_info.txt file, containing names of lanes/chips    
`./00-scripts/00_prepare_lane_info.sh`

Trim for adapters and too short reads    
`./00-scripts/01_cutadapt.sh <numCPUs>`    
Record the results of the trimming for reference

Run fastqc on trimmed data and summarize results     
`mkdir 02-raw/trimmed/fastqc_trimmed/`    
`multiqc -o 02-raw/trimmed/fastqc_trimmed/ 02-raw/trimmed/fastqc_trimmed`         

Prepare a sample information file manually, using the example file, sample_information.csv as a template.    


### De-multiplex your data    
Trim with process_radtags
`./00-scripts/02_process_radtags.sh 70 nsiI mspI`

Rename and copy
`./00-scripts/03_rename_samples.sh`

Determine number of reads per sample
`cd 04-all_samples`    
for i in $(ls *.fq.gz) ; do echo $i ; gunzip -c $i | grep -cE '^@' - ; done
# move to excel and plot

# THIS IS WHERE de novo AND REFERENCE-BASED DIVERGE
#
### de novo ###

# assess output of ustacks
grep -E 'Sample ID|Loaded|merged into|Number of utilized reads|After remainders merged' 10-log_files/2017-05-31_11h50m34s_stacks_1a_ustacks.log > 10-log_files/ustacks_output_trimmed.txt
# this provides:
# sample ID; Number RAD-Tags; Number of loci; Coverage stats; Number used reads
# Then run the following assessment script to collect data
./00-scripts/assessment_scripts/01_assess_ustacks.sh
# Then use the R script to get stats and figures on this output data  
/Users/wayne/Documents/miller/Moore_oyster/04_analysis/stacks_workflow_2017-05-31/assessing_results/assessing_results.R

# skip ahead to cstacks


### reference-based ###
# Align against reference genome
# Two options: Genbank genome (longer) or PAG genome (in chromosomes)

# Locations of genomes:
#PAG
/Users/wayne/Documents/miller/Moore_oyster/00_resources/Cgigas_genome_PAG/Assembly_1/C_gigas_assembly_1_all_chr.fa
# NCBI
/Users/wayne/Documents/z-genome-Cgig/GCA_000297895.1_oyster_v9_genomic.fna.gz

# First try with PAG genome
# first need to edit the alignment script to update the variables GENOMEFOLDER and GENOME
00-scripts/bwa_mem_align_reads.sh 3 
# also see lots of alignment settings and mapping settings for more information on what is being retained
# There are options for the PAG genome or the NCBI genome

### Assessing Mapping Results ###
# Then check for things like how many mappings were successful
# number of mappings for one sample:
samtools view ./Barnes_95.bam | wc -l
# number of unique scaffolds being mapped to:
samtools view Barnes_95.bam | awk '{print $3}' - | sort -n | uniq -c | sort -n | wc -l


### Prepare for Stacks ###
# Finally, prepare the population map file, this will go into 01-info_files/population_map.txt
./00-scripts/04_prepare_population_map.sh


### Running stacks with individual components ###
./00-scripts/stacks_1b_pstacks.sh

# Obtain info 
# grep -E 'Alignments|^Kept' 10-log_files/*_1b_pstacks.log > ./10-log_files/pstacks_output.txt

# cstacks
./00-scripts/stacks_2_cstacks.sh

# Obtain info 
# grep -E 'newly added|loci in the catalog| were matched to a catalog locus\.' 10-log_files/2017-05-31_14h41m03s_stacks_2_cstacks.log > 10-log_files/cstacks_output_trimmed.txt

# run assessment script
#00-scripts/assessment_scripts/02_assess_cstacks.sh
# then go to the same R script to make figures

./00-scripts/stacks_3_sstacks.sh
# see the output by:
# grep -E 'Processing|matching loci' 10-log_files/*sstacks.log | grep -v 'haplotypes examined' - > 10-log_files/sstacks_output_trimmed.txt
# this gives you how many loci were matched against the catalog per sample

# populations (basically just to get a vcf, as filtering done in next step)
./00-scripts/stacks_4_populations.sh

# Filtering
00-scripts/05_filter_vcf.py -i 05-stacks/batch_1.vcf -o 05-stacks/batch_1_filt.vcf -c 1 -m 4 -I 8 -p 70 --use_percent -a 0.01 -A 0.05 -s 20 -H 0.6

# Limit to a single SNP (max_maf)
00-scripts/utility_scripts/extract_snp_with_max_maf.py 05-stacks/batch_1_filt.vcf 05-stacks/batch_1_filt_max_maf.vcf

# How many loci in your vcf file? 
grep -vE '#' 05-stacks/batch_1_filt_max_maf.vcf  | wc -l

# What is the total number of reads present in your output vcf?
vcftools --vcf ./05-stacks/batch_1_filt_max_maf.vcf --site-depth --out 05-stacks/batch_1_filt_max_maf_site-depth
# Then sum up the column in excel

# Run populations again on your filtered vcf (not the single snp one)
populations --in_vcf 05-stacks/1M_filt.vcf --fstats -f p_value --out_path ./05-stacks/re-run_popn_1M/ -M 01-info_files/population_map.txt 
