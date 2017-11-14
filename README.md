# MOP oyster population genetics 
Instruction guide and scripts for the Moore Oyster Project genotyping using Eric Normandeau's stacks_workflow pipeline    
https://github.com/enormandeau/stacks_workflow

Requirements:    
`cutadapt`    
`fastqc`   
`multiqc`   
`bwa`   
`samtools`    
`stacks`    


## Prepare Data

### Setup 

1. Put all raw data in `02-raw` using cp or cp -l    

2. Manually prepare a sample info file (see example file sample_information.csv).   
This must be a tab-delimited text file, but the file name is .csv.    

3. Add the barcodes file that contains all barcodes for the project. This corresponds to the sample information file.    

### Cleanup

Run fastqc and summarize results    
```
mkdir 02-raw/fastqc_raw    
fastqc 02-raw/*.fastq.gz -o 02-raw/fastqc_raw/ -t 5    
multiqc -o 02-raw/fastqc_raw/ 02-raw/fastqc_raw   
```

Use automated script to prep lane_info.txt file    
`./00-scripts/00_prepare_lane_info.sh`

Trim for adapters and too short reads    
`./00-scripts/01_cutadapt.sh <numCPUs>`    
Record the results of the trimming for reference

Run fastqc on trimmed data and summarize results     
```
mkdir 02-raw/trimmed/fastqc_trimmed/    
multiqc -o 02-raw/trimmed/fastqc_trimmed/ 02-raw/trimmed/fastqc_trimmed       
```

### De-multiplex
Trim with process_radtags    
`./00-scripts/02_process_radtags_2_enzymes.sh 80 nsiI mspI` 

Use automated script to rename and copy samples    
`./00-scripts/03_rename_samples.sh`

Prepare for Stacks    
Use automated script to prepare the population map file     
`./00-scripts/04_prepare_population_map.sh`

### Align samples against reference genome
Here use C. gigas assembly from NCBI    
Index the reference genome for use with bwa   
`bwa index GCA_000297895.1_oyster_v9_genomic.fna.gz`

Edit the following script to point to the correct GENOMEFOLDER (full path!) and GENOME variables, then run it        
`00-scripts/bwa_mem_align_reads.sh 6`     

Compare total number of reads to the total number of mappings per sample using the automated script:
`./../ms_oyster_popgen/01_scripts/assess_results.sh`    
This produces files in `04-all_samples`: `reads_per_sample_table.txt` and `mappings_per_sample.txt` and a graph in the main directory of number reads per sample vs number of mappings.

Compare total number of reads per sample to the number of unique scaffolds being mapped against using the script:    
`./../ms_oyster_popgen/01_scripts/determine_number_unique_scaff_mapped.sh`    
This produces a graph as well as some summary statistics.   


### Remove individuals with too few reads     
Note: first requires that the following script was run in prev. step:    
`./../ms_oyster_popgen/01_scripts/assess_results.sh`    

Identify samples with less than 1.5 M reads    
`awk '$2 < 1500000 { print $1 } ' 04-all_samples/reads_per_sample_table.txt > 04-all_samples/samples_to_remove_under1.5M.txt`

Make directory to remove too few read indiv.    
`mkdir 04-all_samples/removed_samples`

Move bam files into the removed_samples directory    
`cd 04-all_samples/ ; for i in $(sed 's/\.fq\.gz/\.bam/g' samples_to_remove_under1.5M.txt ) ; do mv $i removed_samples/ ; done ; cd ..`

Also move fq.gz files into the removed_samples directory    
`cd 04-all_samples/ ;  for i in $(cat samples_to_remove_under1.5M.txt) ; do mv $i removed_samples/ ; done ; cd ..`   

Proceed to Stacks steps ahead [pstacks](#stacks-steps)

If you want to know descriptive stats for only the retained samples, re-run the `assess_results.sh` script.   


## Stacks steps
### pstacks
Adjust the number of threads and launch    
`./00-scripts/stacks_1b_pstacks.sh`

Obtain some info on your pstacks alignment results from the log file:   
`./../ms_oyster_popgen/01_scripts/01_assess_pstacks.sh`   
Produces output `output_pstacks_results.csv` and graph with num reads per sample and average locus coverage per sample.

### cstacks
Edit the following script to use the -g flag (use genomic location) and turn off the -n flag (number mismatches allowed), and incr p (threads)   
`./00-scripts/stacks_2_cstacks.sh`

Assessment: per sample, number of loci matched to the catalog, number of loci added to catalog (starts at second sample).    
`./../ms_oyster_popgen/01_scripts/02_assess_cstacks.sh`     
Produces output `cstacks_output_table.csv`

### sstacks
Edit following script to use -g flag
`./00-scripts/stacks_3_sstacks.sh`     
Log file can be viewed to see how many loci are matched against the catalog. 

### populations
Create .vcf with minimal filtering.    
Edit following script to remove the log-likelihood filter which is not working with alignment based data (--lnl_lim)     
`./00-scripts/stacks_4_populations.sh`

### rxstacks
Edit the following script to remove the log-likelihood filter (--lnl_filter ; --lnl_lim -10)    
`00-scripts/stacks_5b_rxstacks.sh`

### cstacks (round 2)
Edit to add -g flag and remove -n flag    
`00-scripts/stacks_6_cstacks_rx.sh`

### sstacks (round 2) 
Edit to add -g flag    
`00-scripts/stacks_7_sstacks_rx.sh`

### populations (round 2)
Edit to remove lnl_lim
`00-scripts/stacks_8_populations_rx.sh`    
When creating capture panel, also set: min_maf=0.01, p=<the # of pops>, r = 0.4 or 0.5      
When going forward to do pop gen work, set the r value above to 0.7   

### Extra filtering
`00-scripts/05_filter_vcf.py -i 06-stacks_rx/batch_1.vcf -o 06-stacks_rx/batch_1_filt.vcf -c 1 -m 4 -I 8 -p 40 --use_percent -a 0.01 -s 20 -H 0.5 -C 200`     
This provides additional filters of allelic imbalance (-I), max SNPs per locus, max heterozygosity (-H), and max depth (-C)     

### Graph output 
Unfiltered:    
`./00-scripts/05_filter_vcf.py -i 05-stacks/batch_1.vcf -o graphs_before_filters_oyster -g`

Filtered:    
`./00-scripts/05_filter_vcf.py -i 05-stacks/batch_1_filt.vcf -o graphs_after_filters_oyster -g`

Combine:    
`00-scripts/utility_scripts/combine_distribution_graphs.py graphs_before_filters_oyster graphs_after_filters_oyster graphs_both_oyster`

Proceed to:   
[Evaluate number of reads used in output](#evaluate-number-of-reads-used-in-output)    
[RAD Capture panel](#generate-rapture-panel)    
[Final Output](#final-output)    

## Evaluate number of reads used in output    
This will provide information on how many of your reads were actually used in the final output.   
Limit to a single SNP (max_maf) to only count once per locus    
`00-scripts/utility_scripts/extract_snp_with_max_maf.py 06-stacks_rx/batch_1_filt.vcf 06-stacks_rx/batch_1_filt_max_maf.vcf`

Count number of loci in vcf    
`grep -vE '^#' 06-stacks_rx/batch_1_filt_max_maf.vcf  | wc -l`

Count how many reads are being used in your single SNP vcf     
`vcftools --vcf ./06-stacks_rx/batch_1_filt_max_maf.vcf --site-depth --out 06-stacks_rx/batch_1_filt_max_maf_site-depth`

Sum reads     
`grep -vE '^CHROM' 06-stacks_rx/batch_1_filt_max_maf_site-depth.ldepth | awk '{ print $3 }' - | paste -sd+ - | bc`

Count how many reads are in read input file    
`awk '{ print $2 }' 04-all_samples/reads_per_sample_table.txt | paste -sd+ - | bc`

Note (#todo) - make sure this is being done on only the samples that were used after removing samples with too few reads.   

Evaluate total number of mappings    
`awk '{ print $2 }' 04-all_samples/mappings_per_sample_table.txt | paste -sd+ - | bc`


## Final output
This is used to generate fstats on a filtered vcf.   
todo: Needs to be updated    
`populations --in_vcf 06-stacks_rx/1M_filt.vcf --fstats -f p_value --out_path ./06-stacks_rx/re-run_popn_1M/ -M 01-info_files/population_map.txt`

## Generate Rapture panel
This is used to generate a set of tags to be provided for RAD capture bait design.   
Follow these steps:    
* populations module to output batch_1.vcf
* 05-filter_vcf.py at 50% presence to output batch_1_filt.vcf

### Generate single accession output from total fasta with whitelist    
* Obtain a single SNP per locus with max_maf filter    
`00-scripts/utility_scripts/extract_snp_with_max_maf.py 06-stacks_rx/batch_1_filt_p50.vcf 06-stacks_rx/batch_1_filt_p50_max_maf.vcf`
* Identify the catalog locus IDs from the max_maf.vcf    
`grep -vE '^#' 06-stacks_rx/batch_1_filt_p50_max_maf.vcf | awk ' { print $3 } ' - | awk -F_ ' { print $1 } ' - > 06-stacks_rx/whitelist_denovo_max_maf_p50_SNP.txt`
* Edit populations script to use -W flag and point towards the whitelist. Turn on .fasta output and turn off vcf output.    
* Obtain a single Allele 0 record per locus from the fasta output    
`grep -E '^>' 06-stacks_rx/batch_1.fa | awk -FSample_ '{ print $1 }' - | uniq > 06-stacks_rx/obtain_one_record_per_accn_list.txt`
* Use this record list to obtain the single record:    
`while read p; do grep -A1 -m1 $p".*Allele_0" 06-stacks_rx/batch_1.fa ; done < 06-stacks_rx/obtain_one_record_per_accn_list.txt > 06-stacks_rx/batch_1_filtered_single_record.fa`    



## fineRADstructure
Use the output in `06-stacks_rx/batch_1.haplotypes.tsv`    
`cut -f1 --complement batch_1.haplotypes.tsv | cut -f1 --complement - | grep -v 'consensus' | sed 's/-//g' - > batch_1.haplotypes_for_export.tsv`

Need to install `ape` and `XML` 
ape can be installed within R
To install XML, need to first do: 
`sudo apt-get install libxml2-dev`
