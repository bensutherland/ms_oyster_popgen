### Analyze MOP oyster populations using Eric Normandeau's stacks_workflow pipeline
https://github.com/enormandeau/stacks_workflow

Requirements:    
`cutadapt`    
`fastqc`   
`multiqc`   
`bwa`   
`samtools`    
`stacks`    


### Setup
All raw data in `02-raw` using cp or cp -l    

Manually prepare a sample info file (see example file sample_information.csv).   
This must be a tab-delimited text file, but the file name is .csv.    

Add the barcodes file that contains all barcodes for the project. This corresponds to the sample information file.    

### Cleanup

Run fastqc and summarize results    
`mkdir 02-raw/fastqc_raw`    
`fastqc 02-raw/*.fastq.gz -o 02-raw/fastqc_raw/ -t 5`    
`multiqc -o 02-raw/fastqc_raw/ 02-raw/fastqc_raw`   

Use automated script to prep lane_info.txt file
`./00-scripts/00_prepare_lane_info.sh`

Trim for adapters and too short reads    
`./00-scripts/01_cutadapt.sh <numCPUs>`    
Record the results of the trimming for reference

Run fastqc on trimmed data and summarize results     
`mkdir 02-raw/trimmed/fastqc_trimmed/`    
`multiqc -o 02-raw/trimmed/fastqc_trimmed/ 02-raw/trimmed/fastqc_trimmed`       

### De-multiplex
Trim with process_radtags
`./00-scripts/02_process_radtags_2_enzymes.sh 70 nsiI mspI` 

Use automated script to rename and copy samples    
`./00-scripts/03_rename_samples.sh`


## Reference-based stacks
### Align samples against reference genome
Here use C. gigas assembly from NCBI    
Index the reference genome for use with bwa   
`bwa index GCA_000297895.1_oyster_v9_genomic.fna.gz`

Edit the following script to point to the correct GENOMEFOLDER and GENOME variables, then run it        
`00-scripts/bwa_mem_align_reads.sh 6`     

Compare total number of reads to the total number of mappings per sample using the automated script:
`./../ms_oyster_popgen/01_scripts/assess_results.sh`    
This produces files in `04-all_samples`: `reads_per_sample_table.txt` and `mappings_per_sample.txt`.   
Also produces a graph in the main directory of number reads per sample vs number of mappings.

Compare total number of reads per sample to the number of unique scaffolds being mapped against using the script:    
`./../ms_oyster_popgen/01_scripts/determine_number_unique_scaff_mapped.sh`    
This produces a graph as well as some summary statistics.   

#todo: retain summary statistics for both above in a log file     


### Prep for Stacks
Use automated script to prepare the population map file
`./00-scripts/04_prepare_population_map.sh`

### pstacks
`./00-scripts/stacks_1b_pstacks.sh`

Obtain some info on your pstacks alignment results from the log file:   
`./../ms_oyster_popgen/01_scripts/01_assess_pstacks.sh`   
This provides an output file `output_pstacks_results.csv` that provides sample names, the number of loci in the sample, and the average coverage. Also gives graph.


# cstacks
First edit the cstacks script to enable use of -g for using genomic location rather than sequence similarity.
`./00-scripts/stacks_2_cstacks.sh`

Assess results to get the number of loci matched to catalog and number new loci added per sample for all of the samples (starting with the second sample because there can be no matched for sample 1).
`./../ms_oyster_popgen/02_assess_cstacks.sh`


# sstacks
Edit to use more cores and to use -g flag
`./00-scripts/stacks_3_sstacks.sh`
The log file can be viewed to see how many loci are matched against the total number of loci in the catalog (if needed).

# populations
This step is basically just to create the vcf file, as filtering is done later
`./00-scripts/stacks_4_populations.sh`
Make sure to change the log likelihood to turn it off, as this doesn't work with my current data.
--lnl_lim 

# Filtering
00-scripts/05_filter_vcf.py -i 05-stacks/batch_1.vcf -o 05-stacks/batch_1_filt.vcf -c 1 -m 4 -I 8 -p 70 --use_percent -a 0.01 -A 0.05 -s 20 -H 0.6


### To find number of reads used ###
# Limit to a single SNP (max_maf)
00-scripts/utility_scripts/extract_snp_with_max_maf.py 05-stacks/batch_1_filt.vcf 05-stacks/batch_1_filt_max_maf.vcf

# How many loci in your vcf file? 
grep -vE '^#' 05-stacks/batch_1_filt_max_maf.vcf  | wc -l

# What is the total number of reads present in your output vcf?
`vcftools --vcf ./05-stacks/batch_1_filt_max_maf.vcf --site-depth --out 05-stacks/batch_1_filt_max_maf_site-depth`

# Then sum up the column
`grep -vE '^CHROM' 05-stacks/batch_1_filt_max_maf_site-depth.ldepth | awk '{ print $3 }' - | paste -sd+ - | bc`

# Compare to the total number of reads input
`awk '{ print $2 }' 04-all_samples/reads_per_sample_table.txt | paste -sd+ - | bc`

# And the number of mappings
`awk '{ print $2 }' 04-all_samples/mappings_per_sample_table.txt | paste -sd+ - | bc`

# Run populations again on your filtered vcf (not the single snp one)
populations --in_vcf 05-stacks/1M_filt.vcf --fstats -f p_value --out_path ./05-stacks/re-run_popn_1M/ -M 01-info_files/population_map.txt 



### Graphing ###
Before filters
`./00-scripts/05_filter_vcf.py -i 05-stacks/batch_1.vcf -o graphs_before_filters_oyster -g`

After filters
`./00-scripts/05_filter_vcf.py -i 05-stacks/batch_1_filt.vcf -o graphs_after_filters_oyster -g`

Combine
`00-scripts/utility_scripts/combine_distribution_graphs.py graphs_before_filters_oyster graphs_after_filters_oyster graphs_both_oyster`


### Remove low read individuals ###
# identify those with less than 1.5 M reads
`awk '$2 < 1500000 { print $1 } ' 04-all_samples/reads_per_sample_table.txt`

# Save as a file
`awk '$2 < 1500000 { print $1 } ' 04-all_samples/reads_per_sample_table.txt > 04-all_samples/samples_to_remove_under1.5M.txt`

mkdir 04-all_samples/removed_samples

cd 04-all_samples/ ; for i in $(sed 's/\.fq\.gz/\.bam/g' samples_to_re
move_under1.5M.txt ) ; do mv $i removed_samples/ ; done ; cd ..


