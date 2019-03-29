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
`fineRADstructure`    http://cichlid.gurdon.cam.ac.uk/fineRADstructure.html
`ape`   (install within R)   
`libxml2-dev`   (req'd for XML) `sudo apt-get install libxml2-dev`
`XML`    
`fastStructure`    http://rajanil.github.io/fastStructure/

## Prepare Data

### Setup 

1. Put all raw data in `02-raw` using cp or cp -l    

2. Manually prepare a sample info file (see example file sample_information.csv).   
This must be a tab-delimited text file, but the file name is .csv.    

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
fastqc -t 5 02-raw/trimmed/*.fastq.gz -o 02-raw/trimmed/fastqc_trimmed/
multiqc -o 02-raw/trimmed/fastqc_trimmed/ 02-raw/trimmed/fastqc_trimmed       
```

### De-multiplex
Trim with process_radtags in parallel (here 8 cores)     
`00-scripts/02_process_radtags_2_enzymes_parallel.sh 80 nsiI mspI 8`    

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

Want to see how many reads are present in total?    
`awk '{ print $2 } ' 04-all_samples/reads_per_sample_table.txt | paste -sd+ - | bc`

And before assigning (must divide by 4) :   
`for i in $(ls 02-raw/*.fastq.gz) ; do echo $i ; gunzip -c $i | wc -l ; done`

Compare total number of reads per sample to the number of unique scaffolds being mapped against using the script:    
`./../ms_oyster_popgen/01_scripts/determine_number_unique_scaff_mapped.sh`    
This produces a graph as well as some summary statistics.   

### Remove problematic individuals
Make directory to remove problematic or off-project samples.    
`mkdir 04-all_samples/removed_samples`

This individual was found to have very low mapping relative to read counts, so remove it:   
`mv 04-all_samples/PIP_631.* 04-all_samples/removed_samples/`    

### Remove individuals with too few reads     
Note: this step requires that you have already 'assessed results' above.      

Identify samples with less than set number of reads    
`awk '$2 < 1000000 { print $1 } ' 04-all_samples/reads_per_sample_table.txt > 04-all_samples/samples_to_remove.txt`

Move bam files into the removed_samples directory    
`cd 04-all_samples/ ; for i in $(sed 's/\.fq\.gz/\.bam/g' samples_to_remove.txt ) ; do mv $i removed_samples/ ; done ; cd ..`

Also move fq.gz files into the removed_samples directory    
`cd 04-all_samples/ ;  for i in $(cat samples_to_remove.txt) ; do mv $i removed_samples/ ; done ; cd ..`   

Retain the original report files (e.g. reads, mappings):    
`mv 04-all_samples/samples_to_remove.txt 04-all_samples/mappings_per_sample_table.txt 04-all_samples/reads_per_sample_table.txt ./reads_and_mappings_current.pdf 04-all_samples/removed_samples/`

Recalculate sample stats with the limited number of samples:     
`./../ms_oyster_popgen/01_scripts/assess_results.sh`    

### Determine how many samples you have per population
Optional:      
`ls -1 04-all_samples/*.bam | awk -F"/" '{ print $2 }' - | awk -F"_" '{ print $1 }' - | sort -n | uniq -c`    

### Remove filtered individuals from the population map
Use the following to only keep retained samples in your population map:    
`../ms_oyster_popgen/01_scripts/redo_population_map.sh`      

This creates the file `01-info_files/population_map_retained.txt`.     

### Merge closely-related collections in the population map
Optional:      
If you want to merge closely-related populations prior to stacks, use the following:      
`rename_populations_in_map.R`

This will produce `"01-info_files/population_map_retained_renamed.txt"`.    

At the final stage of stacks, when filtering out variants using the populations module, this new file is useful.      

## Stacks steps
Note: we now have a runall option that will run through all stacks steps, using the parameters set in each individual script.
`./../ms_oyster_popgen/01_scripts/runall_stacks.sh`

### pstacks
Adjust the number of threads and launch    
`./00-scripts/stacks_1b_pstacks.sh`

Optional: obtain some info on your pstacks alignment results from the log file:   
`./../ms_oyster_popgen/01_scripts/01_assess_pstacks.sh`   
Requires on an up-to-date assess results calculation.    
This produces `output_pstacks_results.csv` and a graph with num reads per sample and average locus coverage per sample.     

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
`00-scripts/stacks_8_populations_rx.sh`    

Current parameters:     
```
p=<number pops>     
r=0.7    
min_maf=0.01
#lnl_lim # remove if using -g flag
```

Need to know how many populations remain?    
`awk -F_ '{ print $1 }' 01-info_files/population_map_retained.txt | sort -n | uniq -c | less`    

If building RAD capture panel, more loci would be retained with r = 0.4 or 0.5

First set to `write_single_snp`' to export:   
- VCF for use in diversity stats [Diversity](#nucleotide-diversity)     
- plink ped/map files for use in adegenet [General Stats](#hierfstat-and-adegenet)

```
# Save output to analysis folders
mkdir 09-diversity_stats    
mv 06-stacks_rx/batch_1.vcf 09-diversity_stats
mkdir 11-adegenet_analysis
mv 06-stacks_rx/batch_1.ped 06-stacks_rx/batch_1.map 11-adegenet_analysis/
```

Then re-run populations module with multiple SNP allowed and haplotype VCF activated to:    
- haplotypes VCF for fineRADstructure [fineRADstructure](#fineradstructure)    
```
# Save haplotype output to analysis folder
mkdir 08-fineRADstructure
mv 06-stacks_rx/batch_1.haplotypes.vcf  08-fineRADstructure 
```


## hierfstat and adegenet
Use a single SNP per locus for this, and output in plink format, which will allow you to move genomic SNP data into adegenet.    

Generate input file for adegenet:    
`plink --ped 06-stacks_rx/batch_1.plink.ped --map 06-stacks_rx/batch_1.plink.map --maf 0.01 --recodeA --noweb --out 06-stacks_rx/batch_1`

Go to the Rscript `01_scripts/adegenet.R` and load in the .raw file from plink.      
This script will allow you to import the data, build a neighbour-joining tree, perform PCA, DAPC, and bootstrap Fst values for all populations. The script will end with a conversion to a genind object from genlight, and saving out as `11-other_stats/adegenet_output.RData`.    

Next, use the script `01_scripts/adegenet_sep_pops.R` to import the genind object from above, then perform various stats on separated subcomponents of the full dataset. For example:    
* All local samples from BC
* Spatial vs. temporal in BC (Hisnit and Serpentine)
* Selection experiment BC
* Selection experiment France
* Selection experiment China

Result figures and tables will be labeled using their respective subcomponent, and placed in the `11-other_stats` folder.     

## Nucleotide diversity
Use a single SNP per locus for this, and move the output vcf to a new directory.     
```
mkdir 09-diversity_stats
cp 06-stacks_rx/batch_1.vcf 09-diversity_stats/
```

Use automated script to get sample names for each data subset, then calculate genetic diversity (Pi) for each subset with vcftools.       
The script will generate a per locus pi value for the following sets of samples:      
1. Hisnit/Pendrell/Pipestem/Serpentine 
2. BC wild-to-farm (PENF)
3. Rosewall
4. DeepBay
5. FranceW
6. FranceC
7. China Farm 1 (QDC)
8. China Farm 2 (RSC)
9. China Wild (CHN)
10. China wild-to-farm (CHNF)
11. Guernsey
12. Japan        
`../ms_oyster_popgen/01_scripts/nuc_diversity.sh`

Also obtain Fis per individual:    
```
vcftools --vcf 09-diversity_stats/batch_1.vcf --het    
mv out.het 09-diversity_stats/batch_1.het
```

Then use the RScript `diversity_comparison.R`    

## Haplotype Analysis
### fineRADstructure
Use multiple SNP per locus for fineRADstructure input (`06-stacks_rx/batch_1.haplotypes.tsv`)   

Make a new directory, move the haplotypes text file to this new directory, and format the file as an input to fineRADstructure:    
```
mkdir 08-fineRADstructure
cp 06-stacks_rx/batch_1.haplotypes.tsv 08-fineRADstructure
# It is best to use the Stacks2fineRAD.py script to prepare stacks output into fineRADstructure input as this gives some info about the loci and removes instances of triploid alleles.   
python /home/ben/Programs/fineRADstructure/Stacks2fineRAD.py -i 08-fineRADstructure/batch_1.haplotypes.tsv -n 10 -m 30
# This will output a file that has loci filtered based on ploidy, and samples filtered by missing data. 
```

If you want to drop individuals from this analysis to reduce the total size of the graphs, this is best done after the haplotypes.tsv has been translated to the fiiltered fineRADpainter input file (i.e. after the Stacks2fineRAD step). Then one can use `fineRADstructure_input_reduction.R`. PLEASE NOTE: this Rscript will write over the original input of the Stacks2fineRAD step.    

In general, follow instructions from the fineRADstructure tutorial (http://cichlid.gurdon.cam.ac.uk/fineRADstructure.html), but in brief:  
1) Calculate co-ancestry matrix:     
` RADpainter paint 08-fineRADstructure/batch_1.haplotypes.tsv.fineRADpainter.lociFilt.samples30%missFilt.txt`
2) Assign individuals to populations:     
`finestructure -x 100000 -y 100000 08-fineRADstructure/batch_1.haplotypes.tsv.fineRADpainter.lociFilt.samples30%missFilt_chunks.out 08-fineRADstructure/batch_1.haplotypes.tsv.fineRADpainter.lociFilt.samples30%missFilt_chunks.mcmc.xml`     
note: this uses by default a Markov chain Monte Carlo (mcmc) without a tree. By default it assumes the data comes from one population (-I 1), and uses 100000 burn in (-x) and sample iterations (-y) for the MCMC.    
3) Build tree: `finestructure -m T -x 10000 08-fineRADstructure/batch_1.haplotypes.tsv.fineRADpainter.lociFilt.samples30%missFilt_chunks.out 08-fineRADstructure/batch_1.haplotypes.tsv.fineRADpainter.lociFilt.samples30%missFilt_chunks.mcmc.xml 08-fineRADstructure/batch_1.haplotypes.tsv.fineRADpainter.lociFilt.samples30%missFilt_chunks.mcmcTree.xml`     

Then plot using the Rscripts available from the following website shown above:    
`fineRADstructurePlot.R` (follow instructions here)    
`FinestructureLibrary.R`
This will produce plots in the working directory.    


## Other Utilities 
Proceed to:   
[Evaluate positions of SNPs in tags](#evaluate-positions-of-snps-in-tags)    
[Evaluate number of reads used in output](#evaluate-number-of-reads-used-in-output)    
[Obtain a single read per locus](#obtain-a-single-read-per-locus)    

## Evaluate position of SNPs in tags
Curious as to where the SNPs are in your tags?    
`grep -vE '^#' 06-stacks_rx/batch_1_filt.vcf | awk '{ print $3 }' - | awk -F_ '{ print $2 }' - | sort -n > snp_positions_only_positions.txt`    
Then use the script `snp_position_plot.R` to plot a histogram of your data.   


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

Note: make sure to calculate this on only the samples used after removing those with too few reads.   

Evaluate total number of mappings    
`awk '{ print $2 }' 04-all_samples/mappings_per_sample_table.txt | paste -sd+ - | bc`

## Obtain a single read per locus
Use a whitelist to generate a single accession per locus fasta file
* Use the vcf that only has a single SNP per locus    
* Identify the catalog locus IDs from the vcf (3rd column, first unit (second unit is the position of the SNP in the tag))    
`grep -vE '^#' 06-stacks_rx/batch_1.vcf | awk ' { print $3 } ' - | awk -F_ ' { print $1 } ' - > 06-stacks_rx/whitelist.txt`
* Edit populations to use -W flag and point towards the whitelist. Turn on .fasta output and turn off vcf output. (todo: give a one-liner for this instead of using script)    
* Obtain a single Allele 0 record per locus from the fasta output    
`grep -E '^>' 06-stacks_rx/batch_1.fa | awk -FSample_ '{ print $1 }' - | uniq > 06-stacks_rx/obtain_one_record_per_accn_list.txt`
* Obtain the single record:    
`while read p; do grep -A1 -m1 $p".*Allele_0" 06-stacks_rx/batch_1.fa ; done < 06-stacks_rx/obtain_one_record_per_accn_list.txt > 06-stacks_rx/batch_1_filtered_single_record.fa`    

## Identify position of loci in a reference genome
e.g. find loci that are high loading for dapc1
`awk -F"," ' $2 > 0.001 { print $1 }' 11-other_stats/locus_and_dapc1_var_contrib.csv | grep -vE '^mname' | awk -F"_" '{ print $1"_"$2 }' - > 11-other_stats/high_dapc_var_loci_to_keep.csv`

This can be used as a white list to export a fasta and get a single record per locus as described above from the populations module. This will produce `06-stacks_rx/batch_1_filtered_single_record.fa` 
This type of file is also produced by the previous section 'Obtain a single read per locus', and this probably makes more sense to use, as there is no point in doing this multiple times, just do it for all markers at once.    

Align against the reference genome:    
`bwa mem -t 6 ~/Documents/genomes/C_gigas_assembly_1_all_chr.fa 06-stacks_rx/batch_1_filtered_single_record.fa > 06-stacks_rx/batch_1_filtered_single_record.sam`    
`samtools view -Sb -F 4 -q 5 06-stacks_rx/batch_1_filtered_single_record.sam > 06-stacks_rx/batch_1_filtered_single_record_q.bam`

Obtain info from the alignment
`samtools view 06-stacks_rx/batch_1_filtered_single_record_q.bam | awk '{ print $1 "," $3 "," $4 }' - > 11-other_stats/aligned_marker_info.csv`


## Run fastStructure
Run plink to convert plink output from stacks into a .bed, .bim, and .fam file for fastStructure input     
`plink --noweb --file 06-stacks_rx_output_296_samples_many_sep_pops/batch_1.plink --make-bed --out 06-stacks_rx_output_296_samples_many_sep_pops/batch_1`

Run faststructure:    
`python ~/Programs/fastStructure/structure.py -K 5 --input=06-stacks_rx_output_296_samples_many_sep_pops/batch_1 --output=test_structure/`

Get labels for distruct:    
`awk '{ print $2 }' 06-stacks_rx_output_296_samples_many_sep_pops/batch_1.fam > structure_labels.txt`     

Run distruct.py to plot results:     
`python ~/Programs/fastStructure/distruct.py -K 5 --input=test_structure --output=test_structure_out.svg --popfile=structure_labels.txt`

Test multiple values of K:    
`for i in $(seq 1 10) ; do echo $i ; python ~/Programs/fastStructure/structure.py -K $i --input=06-stacks_rx_output_296_samples_many_sep_pops/batch_1 --output=test_structure/test_structure ; done`

Evaluate:    
`python ~/Programs/fastStructure/chooseK.py --input=test_structure/test_structure`
