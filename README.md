# MOP oyster population genetics 
Instruction guide and scripts using Eric Normandeau's stacks_workflow pipeline    
https://github.com/enormandeau/stacks_workflow

Requirements:    
`cutadapt`    
`fastqc`   
`multiqc`   
`bwa`   
`samtools`    
`stacks`    


## Setup
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

Prepare for Stacks
Use automated script to prepare the population map file
`./00-scripts/04_prepare_population_map.sh`


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

## Stacks steps
### pstacks
`./00-scripts/stacks_1b_pstacks.sh`
For unmapped reads, 

Obtain some info on your pstacks alignment results from the log file:   
`./../ms_oyster_popgen/01_scripts/01_assess_pstacks.sh`   
Produces output `output_pstacks_results.csv` and graph with num reads per sample and average locus coverage per sample.

### cstacks
For mapped reads, edit following script to use -g to use genomic location instead of seq similarity.    
For unmapped reads, edit the following script to not use the -g option, and to only use the `*unmapped.bam` files.    
`./00-scripts/stacks_2_cstacks.sh`

For mapped reads, assess results to determine per sample the number of loci matched to the catalog and the number of new loci added. Note: this starts at sample 2.    
`./../ms_oyster_popgen/02_assess_cstacks.sh`

### sstacks
Edit following script to use -g flag and use more cores     
For unmapped reads, edit to not use -g and to only use `*unmapped.bam` files.    
`./00-scripts/stacks_3_sstacks.sh`     
Log file can be viewed to see how many loci are matched against the catalog. 

### populations
Basically only to create a .vcf with minimal filtering. Edit script to remove the log-likelihood filter which is not working with alignment based data (--lnl_lim)     
`./00-scripts/stacks_4_populations.sh`


## Filtering
Use vcf filtering script    
`00-scripts/05_filter_vcf.py -i 05-stacks/batch_1.vcf -o 05-stacks/batch_1_filt.vcf -c 1 -m 4 -I 8 -p 70 --use_percent -a 0.01 -A 0.05 -s 20 -H 0.6`

### Graph output 
Unfiltered:    
`./00-scripts/05_filter_vcf.py -i 05-stacks/batch_1.vcf -o graphs_before_filters_oyster -g`

Filtered:    
`./00-scripts/05_filter_vcf.py -i 05-stacks/batch_1_filt.vcf -o graphs_after_filters_oyster -g`

Combine:    
`00-scripts/utility_scripts/combine_distribution_graphs.py graphs_before_filters_oyster graphs_after_filters_oyster graphs_both_oyster`


## Evaluate number of reads used in output    
Limit to a single SNP (max_maf) to only count once per locus    
`00-scripts/utility_scripts/extract_snp_with_max_maf.py 05-stacks/batch_1_filt.vcf 05-stacks/batch_1_filt_max_maf.vcf`

Count number of loci in vcf    
`grep -vE '^#' 05-stacks/batch_1_filt_max_maf.vcf  | wc -l`

Count how many reads are being used in your single SNP vcf     
`vcftools --vcf ./05-stacks/batch_1_filt_max_maf.vcf --site-depth --out 05-stacks/batch_1_filt_max_maf_site-depth`

Sum reads     
`grep -vE '^CHROM' 05-stacks/batch_1_filt_max_maf_site-depth.ldepth | awk '{ print $3 }' - | paste -sd+ - | bc`

Count how many reads are in read input file    
`awk '{ print $2 }' 04-all_samples/reads_per_sample_table.txt | paste -sd+ - | bc`

Evaluate total number of mappings    
`awk '{ print $2 }' 04-all_samples/mappings_per_sample_table.txt | paste -sd+ - | bc`


## Remove individuals with too few reads     
Identify samples with less than 1.5 M reads    
`awk '$2 < 1500000 { print $1 } ' 04-all_samples/reads_per_sample_table.txt > 04-all_samples/samples_to_remove_under1.5M.txt`

Make directory to remove too few read indiv.    
`mkdir 04-all_samples/removed_samples`

Move bam files into the removed_samples directory    
`cd 04-all_samples/ ; for i in $(sed 's/\.fq\.gz/\.bam/g' samples_to_remove_under1.5M.txt ) ; do mv $i removed_samples/ ; done ; cd ..`

Go back and rerun the Stacks section starting at [pstacks](#stacks-steps)
Once you have run it again, use your final, filtered vcf file in the next stage.

## Final output
`populations --in_vcf 05-stacks/1M_filt.vcf --fstats -f p_value --out_path ./05-stacks/re-run_popn_1M/ -M 01-info_files/population_map.txt`


## Generate Rapture panel by combining de novo and reference alignment
**de novo**    
Follow these steps:    
* populations module to output batch_1.vcf
* filter_vcf.py at 50% presence to output batch_1_filt.vcf
* obtain a single SNP per locus with max_maf filter `00-scripts/utility_scripts/extract_snp_with_max_maf.py 05-stacks/batch_1_filt_p50.vcf 05-stacks/batch_1_filt_p50_max_maf.vcf`
* obtain the catalog locus IDs from the vcf
`grep -vE '^#' 05-stacks/batch_1_filt_p50_max_maf.vcf | awk ' { print $3 } ' - | awk -F_ ' { print $1 } ' - > 05-stacks/whitelist_denovo_max_maf_p50_SNP.txt`
* Edit populations script to use -W flag and point towards the whitelist. Turn on .fasta output.    
* Obtain a single Allele 0 record per locus from the fasta output
`grep -E '^>' 05-stacks/batch_1.fa | awk -FSample_ '{ print $1 }' - | uniq > 05-stacks/obtain_one_record_per_accn_list.txt`
* Use this record list to obtain the single record:    
`while read p; do grep -A1 -m1 $p".*Allele_0" 05-stacks/batch_1.fa ; done < 05-stacks/obtain_one_record_per_accn_list.txt > 05-stacks/batch_1_filtered_single_record.fa`    

**reference-based**
Perform all of the above _de novo_ steps but on the output from the reference-based batch_1.vcf. Use a separate 05-stacks folder.    

### Compare the results
Index the aligned fasta output:    
`bowtie2-build -f batch_1_filtered_single_record_aligned.fa --threads 6 batch_1_filtered_single_record_aligned`

Map de novo against the aligned version    
`bowtie2 -x batch_1_filtered_single_record_aligned -f batch_1_filtered_single_record_denovo.fa --end-to-end --threads 6 > denovo_vs_aligned.sam`    

Obtain unmapped reads from this sam file     
`samtools view -Sf 4 denovo_vs_aligned.sam > denovo_vs_aligned_unmapped.sam`    
`awk '{ print $1 }' denovo_vs_aligned_unmapped.sam > new_markers_to_add.txt`

Extract _de novo_-only loci from the _de novo_ fasta file to create a new fasta file.     
`while read p; do grep -A1 $p ./batch_1_filtered_single_record_denovo.fa ; done < new_markers_to_add.txt > new_markers_to_add.fa`

Concatenate to the genome fasta file, index this for use with bwa, then redo the entire alignment and stacks pipeline with this as the new reference genome.

Issues: 
* too high coverage? Incorporate a screen filter for -C flag in 05-filter_vcf.py


