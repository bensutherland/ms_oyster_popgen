# Moore Oyster Project - Population Genetics
Manuscript analysis instruction guide and associated scripts for the Moore Oyster Project Population Genetics analysis. This pipeline is still in development stage, and comes with no guarantees.            
Primarily uses the repo https://github.com/bensutherland/stacks_workflow, which is a fork of the original repo by Eric Normandeau in Labo Bernatchez https://github.com/enormandeau/stacks_workflow but customized for Stacks v2.0.       


### Requirements    
`cutadapt` https://cutadapt.readthedocs.io/en/stable/index.html    
`fastqc` https://www.bioinformatics.babraham.ac.uk/projects/fastqc/   
`multiqc` https://multiqc.info/   
`bwa` http://bio-bwa.sourceforge.net/   
`samtools (v.1.9)` http://www.htslib.org/    
`stacks (v2.3e)` http://catchenlab.life.illinois.edu/stacks/     
`fineRADstructure` http://cichlid.gurdon.cam.ac.uk/fineRADstructure.html     
`ape` (install within R)     
`libxml2-dev` (req'd for XML, install Linux w/ apt-get)       
`XML`    
`fastStructure`    http://rajanil.github.io/fastStructure/     
`stacks_workflow` original from E. Normandeau, but for Stacks v2 use fork: https://github.com/bensutherland/stacks_workflow        
`pcadapt`       
`bayescan`   cmpg.unibe.ch/software/BayeScan/download.html

All analysis is done within the `stacks_workflow` repo, which should be contained within the same parent directory as this repo (at the same level).       

## 1. Preparing Data
### a. Set up 
1. Put all raw data in `02-raw` using cp or cp -l    
2. Prepare the sample info file (see template in repo sample_information.csv). Note: tab-delimited, even though name is .csv.    
3. Download reference genome (oyster_v9): https://www.ncbi.nlm.nih.gov/assembly/GCF_000297895.1/      
note: source citation is Zhang et al. 2012, Nature. https://www.ncbi.nlm.nih.gov/pubmed/22992520/       


### b. Clean data
View raw data with fastqc and multiqc:    
```
mkdir 02-raw/fastqc_raw    
fastqc 02-raw/*.fastq.gz -o 02-raw/fastqc_raw/ -t 5    
multiqc -o 02-raw/fastqc_raw/ 02-raw/fastqc_raw   
```
Prepare lane_info.txt file with automated script:    
`./00-scripts/00_prepare_lane_info.sh`    

Trim adapters and too short reads:    
`./00-scripts/01_cutadapt.sh <numCPUs>`    

View trimmed data with fastqc and multiqc:     
```
mkdir 02-raw/trimmed/fastqc_trimmed/    
fastqc -t 5 02-raw/trimmed/*.fastq.gz -o 02-raw/trimmed/fastqc_trimmed/
multiqc -o 02-raw/trimmed/fastqc_trimmed/ 02-raw/trimmed/fastqc_trimmed       
```

### c. De-multiplex reads
Detect cut sites and barcodes to de-multiplex and truncate reads to 80 bp with process_radtags in parallel:     
`00-scripts/02_process_radtags_2_enzymes_parallel.sh 80 nsiI mspI 8`    

Collect samples from multiple runs together and rename samples:    
`./00-scripts/03_rename_samples.sh`

Prepare the population map file:     
`./00-scripts/04_prepare_population_map.sh`

## 2. Align samples and quality filter on samples
### a. Align samples against the reference genome
Index the reference genome with bwa, here using _Crassostrea gigas_:    
`bwa index GCA_000297895.1_oyster_v9_genomic.fna.gz`    

Point the alignment script to the reference genome (full path) using the GENOMEFOLDER and GENOME variables, and run:         
`00-scripts/bwa_mem_align_reads.sh 6`     

### b. Inspect alignment results
Compare per-sample reads and alignments, and per-sample reads and number of aligned scaffolds:      
`./../ms_oyster_popgen/01_scripts/assess_results.sh`    
`./../ms_oyster_popgen/01_scripts/determine_number_unique_scaff_mapped.sh`    

Produces:      
```
04-all_samples/reads_per_sample_table.txt
04-all_samples/mappings_per_sample.txt
# A graph of number reads and aligned reads
# a graph of number reads and scaffolds mapped
```

Other calculations:
```
# Total reads in all samples:     
awk '{ print $2 } ' 04-all_samples/reads_per_sample_table.txt | paste -sd+ - | bc
# Total reads before de-multiplexing (note: divide by 4 due to fastqc):   
for i in $(ls 02-raw/*.fastq.gz) ; do echo $i ; gunzip -c $i | wc -l ; done
```

### c. Filter out low reads/alignments individuals
Make directory to remove samples:    
`mkdir 04-all_samples/removed_samples`

1. Remove specific individuals    
```
# Remove an individual with low alignments:   
mv 04-all_samples/PIP_631.* 04-all_samples/removed_samples/   

# Remove duplicate individuals (identified previously through shared microhaps)
mv 04-all_samples/GUR_539.* 04-all_samples/removed_samples/
mv 04-all_samples/PENF_214.* 04-all_samples/removed_samples/

```

2. Remove individuals with too few reads:     
(_note: requires 'Inspect alignment results' output_).      
```
# Identify samples with too few reads (here < 1000000)
awk '$2 < 1000000 { print $1 } ' 04-all_samples/reads_per_sample_table.txt > 04-all_samples/samples_to_remove.txt

# Remove .bam, .fq.gz, and report files from analysis
cd 04-all_samples/ ; for i in $(sed 's/\.fq\.gz/\.bam/g' samples_to_remove.txt ) ; do mv $i removed_samples/ ; done ; cd ..
cd 04-all_samples/ ;  for i in $(cat samples_to_remove.txt) ; do mv $i removed_samples/ ; done ; cd ..
mv 04-all_samples/samples_to_remove.txt 04-all_samples/mappings_per_sample_table.txt 04-all_samples/reads_per_sample_table.txt ./reads_and_mappings_current.pdf 04-all_samples/removed_samples/

# Recalculate sample stats with retained samples
./../ms_oyster_popgen/01_scripts/assess_results.sh
```

### d. Characterize input samples and populations 
```
# Essential steps:    
# Calculate the number of samples per population
ls -1 04-all_samples/*.bam | awk -F"/" '{ print $2 }' - | awk -F"_" '{ print $1 }' - | sort -n | uniq -c > 01-info_files/samples_per_pop.txt      

# Reduce pop map to only include the retained individuals
./../ms_oyster_popgen/01_scripts/redo_population_map.sh
# Produces: 01-info_files/population_map_retained.txt

# Rename pop map to have names instead of numeric codes
`R --slave -e 'source("./../ms_oyster_popgen/01_scripts/translate_popmap_numeric_to_names.R")'`

```


No longer used (originally to merge specif. pops: `./../ms_oyster_popgen/01_scripts/z-unused/rename_populations_in_map.R`)         


## 3. Genotype
### a. Run Stacks v2.0 genotyper:     
Genotype using alignments, with 8 cores:      
`./00-scripts/stacks2_1_gstacks.sh 8`      

### b. Filters and different outputs: 
Create output directories for your data:     
`mkdir 09-diversity_stats 11-adegenet_analysis 08-fineRADstructure`      

To identify loci out of HWE then run multiple outputs, use the following runall script:       
`./../ms_oyster_popgen/01_scripts/runall_outputs.sh`        

This will:     
* Identify loci out of HWE using stacks calculated HWE p-values
* Output single SNP per locus (filtered)
* Output multiple SNP per locus (filtered)

Provides output for:
Diversity analysis: [Diversity Analysis](#nucleotide-diversity)     
in `09-diversity_stats`    

Population differentiation analysis:[General Stats](#hierfstat-and-adegenet)     
in `11-adegenet_analysis`       

Relatedness via shared coancestry with SNPs (related) and haplotypes (fineRADstructure) [fineRADstructure](#fineradstructure)      
in `08-fineRADstructure`


### c. Obtain summary info 
#### Sample size per type    
. # 2019-07-22 todo: automate using Rscript, and output reports not just to screen
`01_scripts/characterize_samples.R`    


## 4. Diversity analysis
### A. Nucleotide diversity
Uses: single SNP VCF        

#### i. Per-site nucleotide diversity
Calculate per-site nucleotide diversity (Pi) for all sites in VCF data subsets with vcftools:       
`./../ms_oyster_popgen/01_scripts/nuc_diversity.sh`     

Data subsets:     
1. BC naturalized: Hisnit/Pendrell/Pipestem/Serpentine 
2. BC naturalized-to-farm (PENF)
3. BC farm (Rosewall)
4. BC farm (DeepBay)
5. France naturalized
6. France naturalized-to-farm (FRAF)
7. China farm 1 (QDC)
8. China farm 2 (RSC)
9. China wild (CHN)
10. China wild-to-farm (CHNF)
11. UK hatchery (Guernsey)
12. Japan wild        

#### ii. Per-individual heterozygosity
Calculate the inbreeding coefficient 'F' for each individual with vcftools:      
`../ms_oyster_popgen/01_scripts/heterozygosity.sh`      

#### iii. Plot results for diversity and heterozygosity
Results will be in `09_diversity_stats`
`R --slave -e 'source("./../ms_oyster_popgen/01_scripts/diversity_comparison.R")'`


## 5. Genetic differentiation and relatedness
### A. Genetic differentiation
Uses: single SNP plink       
Convert from plink to adegenet:     
`plink --ped 11-adegenet_analysis/populations.plink.ped --map 11-adegenet_analysis/populations.plink.map --maf 0.01 --recodeA --noweb --out 11-adegenet_analysis/populations_single_snp_HWE`      
Run adegenet analysis (start):     
`R --slave -e 'source("./../ms_oyster_popgen/01_scripts/adegenet.R")'`       
...which will build a NJ tree, conduct PCA and dAPC, and calculate bootstrapped Fst vals per population pair. Will also prep for data subset-specific analysis         

Run adegenet analysis (data subset-specific):     
`R --slave -e 'source("./../ms_oyster_popgen/01_scripts/adegenet_sep_pops.R")'`       
...which will analyzes the following contrasts:        
* Global analysis (only one pop from BC) `global`    
* BC naturalized analysis `bc_all`    
* BC naturalized spatial vs. temporal (Hisnit and Serpentine) `bc_size`     
* China analysis `ch_all`     
* Parallel naturalized-to-farm (i.e. BC, France, and China) `bc_sel`, `fr_sel`, `chr_sel`
* Domestication (i.e. `DPB vs PEN`, `ROS vs PIP`, `QDC vs CHN`, `RSC vs CHN`, and `GUR vs FRA`    
Results will be in `11_adegenet_analysis`      

Private alleles are identified using the following:     
`../ms_oyster_popgen/01_scripts/private_alleles.R`       
This uses the output from the adegenet analysis above, and will output a graph of private alleles and their individual counts for both a bare minimum analysis and an all-BC compared to farms analysis. (#todo finalize code)     

### B. Relatedness
Note: must run differentiation analysis first.    
#### i. Shared coancestry by related
Input: single SNP, HWE filtered from adegenet output (adegenet_output.RData).      
Open Rscript `01_scripts/relatedeness.R` to translate the genlight obj to related format, and calculate coancestry per population, plotting per-population relatedness metrics.       
Output: See `09-diversity_stats/relatedness_*.pdf`   

#### ii. Shared haplotypes by fineRADstructure
Input: haplotype output from stacks populations as input to fineRADstructure         
Depends that you have run the runall script, and at the haplotype output of the populations module, it will output into SimpleMatrix format.     

In general, follow instructions from the fineRADstructure tutorial (http://cichlid.gurdon.cam.ac.uk/fineRADstructure.html), but in brief:  
1) Calculate co-ancestry matrix:     
`RADpainter paint 21-haplotype_results/populations.haps.radpainter`      
2) Assign individuals to populations:     
`finestructure -x 100000 -y 100000 21-haplotype_results/populations.haps_chunks.out 21-haplotype_results/populations.haps_chunks.out.mcmc.xml`
Uses Markov chain Monte Carlo (mcmc) without a tree, assumes data is from one pop (-I 1).    
3) Build tree:        
`finestructure -m T -x 10000 21-haplotype_results/populations.haps_chunks.out 21-haplotype_results/populations.haps_chunks.out.mcmc.xml 21-haplotype_results/populations.haps_chunks.out.mcmcTree.xml`      

Then plot using the Rscripts adapted from the fineRADstructure site (see above)   
`01_scripts/fineRADstructurePlot.R` (follow instructions here)      
`01_scripts/FinestructureLibrary.R`     
This will produce plots in the working directory.    

## 6. Selection analysis
Input: copy the file populations.snps.vcf to a new folder entitled `13_selection`      
Also make a folder for subsetting samples: `13_selection/sets/`

Prepare a file, within `../ms_oyster_popgen/00_archive/selection_contrasts.txt` that is a space or tab separated text file that shows the contrasts of interest, e.g. 
PEN PENF
QDC CHN

Run the following:     
```
# Identify what samples are in your vcf
bcftools query -l 13_selection/populations.snps.vcf > 13_selection/sets/all_samples.txt

# Create set lists of samples within contrasts of interest using the selection_contrasts.txt file
./../ms_oyster_popgen/01_scripts/set_selection_comparisons.sh
# outputs will be VCF files in 13_selection
```

Use the Rscript to calculate outliers per contrast using pcadapt:     
`01_scripts/outlier_detection_pcadapt.r` 


For bayescan, assuming you have already built your contrast strata files for pcadapt, use the following:     
`./../ms_oyster_popgen/01_scripts/run_bayescan.sh`       
...to make a pop map file, format the bayescan input from the VCF, put the bayescan input into its own folder, then run BayeScan within the folder to produce the standard bayescan output files.    
Output file of note: `<CONTRAST>_fst.txt`


### RDA ###
First need to make a VCF with only the farms:     
Make a directory for your work and change into this directory:     
```
mkdir 13_selection/rda_analysis
cd 13_selection/rda_analysis
cat ../sets/CHN_CHNF.txt ../sets/PEN_PENF.txt ../sets/FRA_FRAF.txt > FARM.txt

# Also do for domestication, but make sure to only keep a single copy of each sample here:     
cat ../sets/DPB_PEN.txt ../sets/GUR_FRA.txt ../sets/QDC_CHN.txt ../sets/ROS_PIP.txt ../sets/RSC_CHN.txt | sort | uniq > DOMESTICATION.txt


# Now select only the FARM.txt samples
vcftools --vcf ./../populations.snps.vcf --keep ./FARM.txt --out FARM --recode
vcftools --vcf ./../populations.snps.vcf --keep ./DOMESTICATION.txt --out DOMESTICATION --recode

# convert to a genotype matrix, one row per individual
vcftools --vcf FARM.recode.vcf --012 --out FARM_geno
vcftools --vcf DOMESTICATION.recode.vcf --012 --out DOMESTICATION_geno

``` 

Then go to `../ms_oyster_popgen/01_scripts/outlier_detection_rda.R`       

This will output plots into a folder `13_selection/rda_analysis`, as well as .csv files that can be used to integrate with the genome plot below. However, for the genome plot you will also need this answer key:     
`grep -vE '^#' 13_selection/populations.snps.vcf | awk '{ print $1 ",", $2 "," $3 }' - | awk -F":" '{ print $1 }' - > 13_selection/chrom_pos_id.csv`    

## 7. Plotting outliers along the chromosome-level assembly
### A. Align loci against the reference genome
The populations run above (single snp) will output a fasta file with a consensus sequence per locus, entitled `09-diversity_stats/populations.loci.fa`. Align this against the most complete reference genome:      

Make a folder to work within:    
`mkdir 12_genome_plot`      

Move marker file to the folder:    
`cp 09-diversity_stats/populations.loci.fa 12_genome_plot`      

Do a bit of an update on the fasta file to make sure that all of your records are unique:    
`cat 12_genome_plot/populations.loci.fa | awk '{ print $1"-"$2 }' - | sed 's/\[//g' | sed 's/\,//g' | sed 's/-$//g' | grep -vE '^#' - > 12_genome_plot/populations.loci_fullname.fa`     

The fasta has all loci, including those that were not passing filters, so use the VCF to identify the names of the filtered loci:      
`cp 09-diversity_stats/populations.snps.vcf 12_genome_plot/`    

Then format the names to match the fullname fasta above:    
`grep -vE '^#' 12_genome_plot/populations.snps.vcf | awk '{ print ">CLocus_" $3 "-" $1 }' -  | awk -F":" '{ print $1 $3}' - | sed 's/\+/\-/g' | sed 's/\-\-/\-/g' > 12_genome_plot/markers_in_vcf.txt`

And use xargs to extract only the relevant lines from the renamed fasta:    
`cat 12_genome_plot/markers_in_vcf.txt | xargs -I{} grep -A1 {} 12_genome_plot/populations.loci_fullname.fa | grep -vE '^--$' - > 12_genome_plot/markers_in_vcf.fa`

Your fasta is now ready for alignment.     


#### Genomic coverage? #### 
How to calculate what proportion of the genome is covered by your markers?    
First, find out the number of bp in your markers present in your filtered VCF:     
`grep -vE '^>' 12_genome_plot/markers_in_vcf.fa | wc -c`      
For example, this is calculated as 1532384 bp (1.5 M bp).     
The oyster genome is assumed to be 600 Mbp.     
Therefore, you have `1532384 / 600000000 * 100` , or 0.255% of your genome represented here.    
#### Alignment ####

First download the assembly from http://dx.doi.org/10.12770/dbf64e8d-45dd-437f-b734-00b77606430a 
Citation: Gagnaire _et al._, 2018. Analysis of Genome-Wide Differentiation between Native and Introduced Populations of the Cupped Oysters _Crassostrea gigas_ and _Crassostrea angulata_, Genome Biology and Evolution, *10*(9), pp. 2518–2534, https://doi.org/10.1093/gbe/evy194 .     

Index the genome:    
`bwa index C_gigas_Assembly_1.fa.gz`      

Align all loci (single consensus locus per marker) against the reference genome:        
```
GENOME="/home/ben/Documents/genomes/C_gigas_Assembly_1.fa.gz" ; bwa mem -t 6 $GENOME 12_genome_plot/markers_in_vcf.fa > 12_genome_plot/markers_in_vcf.sam

# Extract only the markers that map:   
samtools view -Sb -F 4 -q 1 12_genome_plot/markers_in_vcf.sam > 12_genome_plot/markers_in_vcf_filtered.bam

# Obtain alignment info per marker for plotting, saving only the top mapq:   
samtools view 12_genome_plot/markers_in_vcf_filtered.bam | awk '{ print $1 "," $3 "," $4 }' - > 12_genome_plot/aligned_marker_info.csv

# Get a quick overview of how many markers on which chromosomes:    
awk -F, '{ print $2 }' 12_genome_plot/aligned_marker_info.csv | sort | uniq -c

# Sum up the number of alignments
awk -F, '{ print $2 }' 12_genome_plot/aligned_marker_info.csv | sort | uniq -c | awk '{ print $1 }' - | paste -sd+ - | bc

# Compare to the markers going into the alignment:   
wc -l 12_genome_plot/markers_in_vcf.txt

```

The next steps involve:   
1) `../ms_oyster_popgen/01_scripts/gwas_01_load_constant_info.r`     
2) `../ms_oyster_popgen/01_scripts/gwas_02_load_set_contrasts_info.r`     
3) `../ms_oyster_popgen/01_scripts/gwas_03_plot.r`       

This will output a bunch of selection plots in `13_selection`


## Find genomic positions and related genes for top outliers
Here we define top outliers as those that were identified using both RDA and pcadapt. To find those, after all selection scripts have been run, use the following:     
```
# Identify the marker names for common outliers that will correspond to the fasta output 
# FASTA output: 12_genome_plot/markers_in_vcf.fa

# Identify the top outliers: 
cat 13_selection/outliers_common_pcadapt_rda_* | awk '{ print $2 }' - | grep -vE '^matching' - | sort -n | uniq -c | sort -nk1 | awk '{ print ">"$2 }' - > 13_selection/top_outliers.txt

# Use top outliers to extract from FASTA
cat 13_selection/top_outliers.txt | xargs -I{} grep -A1 {} 12_genome_plot/markers_in_vcf.fa > 13_selection/top_outliers.fa

# Go to NCBI page for Pacific oyster genome and manual BLAST using that file to see where these are in the genome. 

```



## Other Utilities 
[Evaluate positions of SNPs in tags](#evaluate-positions-of-snps-in-tags)    
[Evaluate number of reads used in output](#evaluate-number-of-reads-used-in-output)    


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


