# Moore Oyster Project - Population Genetics
Manuscript analysis instruction guide and associated scripts for the Moore Oyster Project Population Genetics analysis.      
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

1. Remove specific individuals (e.g. one with low alignments):   
`mv 04-all_samples/PIP_631.* 04-all_samples/removed_samples/`   

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
# Calculate the number of samples per population
ls -1 04-all_samples/*.bam | awk -F"/" '{ print $2 }' - | awk -F"_" '{ print $1 }' - | sort -n | uniq -c
# Reconstruct population map with only retained individuals 
../ms_oyster_popgen/01_scripts/redo_population_map.sh
# Produces: 01-info_files/population_map_retained.txt
```
Optional to merge closely-related collections in the population map based on previous analyses:     
_Update: this may be no longer necessary because of the two tiers allowed in the map file_
```
./rename_populations_in_map.R
# Produces: 01-info_files/population_map_retained_renamed.txt
# This can be useful for the populations module (below)
```

## 3. Genotype
### a. Run Stacks v2.0 individual steps:     
```
# Genotype using alignments 
./00-scripts/stacks2_1_gstacks.sh
# Filter identified variants
./00-scripts/stacks2_2_populations.sh
```

### Obtain data for different analyses: 
The following steps will provide the outputs needed for downstream analyses. Each individual numbered step will provide a different output needed. All edits required are to successive runs of: `./00-scripts/stacks2_2_populations.sh`       

For simplicity, an output script has been created that will successively run populations with the following changes identified below. Therefore, instead of running the following steps, use the script:      
`./../ms_oyster_popgen/01_scripts/runall_outputs.sh`     

**Make some output directories:**     
```
mkdir 08-fineRADstructure
mkdir 09-diversity_stats
mkdir 11-adegenet_analysis
```

**1. Identify loci out of HWE**
```
# Obtain list of loci out of HWE to use as blacklist
# Toggle ON flag --hwe
# Re-run
./00-scripts/stacks2_2_populations.sh
```

**2. Filter and output single SNP per locus**
```
# Toggle ON blacklist flag (-B) pointing to list of loci out of HWE
# Toggle on single snp flag (--write-single-snp)
# Re-run
./00-scripts/stacks2_2_populations.sh

# Observe number of samples per populations remaining
awk -F_ '{ print $1 }' 01-info_files/population_map_retained.txt | sort -n | uniq -c | less

# Save outputs to analysis folder
cp 07-filtered_vcfs/populations.snps.vcf 09-diversity_stats
cp 07-filtered_vcfs/*plink* 11-adegenet_analysis
```

Provides output for:
Diversity analysis: [Diversity Analysis](#nucleotide-diversity)     
Population differentiation analysis:[General Stats](#hierfstat-and-adegenet)

**3. Filter and output microhaplotypes per locus**
```
# Toggle ON blacklist flag (-B) pointing to list of loci out of HWE
# Toggle OFF single snp flag (--write-single-snp)
# Toggle ON radpainter output flag (--radpainter)

# Compare number SNP vs number loci
grep -vE '^#' 06-stacks_rx/batch_1.vcf | wc -l
grep -vE '^#' 06-stacks_rx/batch_1.haplotypes.vcf | wc -l

# Save outputs to analysis folder
cp 07-filtered_vcfs/populations.haps.vcf 08-fineRADstructure
```
Provides output for:    
Relatedness via shared coancestry with SNPs (related) and haplotypes (fineRADstructure) [fineRADstructure](#fineradstructure)      

## 4. Diversity analysis
### A. Nucleotide diversity
Uses the single SNP VCF moved to `09-diversity_stats` (above).    

#### i. Per-site nucleotide diversity
Use automated script to get sample names for each data subset, then use vcftools to calculate per-site nucleotide diversity (Pi) for all sites in the VCF file for each subset.       
`../ms_oyster_popgen/01_scripts/nuc_diversity.sh`     

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

#### ii. Per-individual heterozygosity
Calculate the inbreeding coefficient 'F' for each individual with vcftools:      
`../ms_oyster_popgen/01_scripts/heterozygosity.sh`      

Then use the RScript `../ms_oyster_popgen/diversity_comparison.R`    


## 5. Genetic differentiation and relatedness
### A. Genetic differentiation
Input: single SNP populations output (HWE filter) in plink format to input to adegenet.    

Single SNP, HWE filter; use plink to convert to adegenet:    
`plink --ped 11-adegenet_analysis/populations.plink.ped --map 11-adegenet_analysis/populations.plink.map --maf 0.01 --recodeA --noweb --out 11-adegenet_analysis/populations_single_snp_HWE`      

Open Rscript `01_scripts/adegenet.R` to begin popgen analysis, loading `11-adegenet_analysis/populations_single_snp_HWE.raw`      
Builds NJ tree, PCA, dAPC, bootstrapped Fst vals per population pair       
Outputs multiple formats (e.g. genlight) as `11-other_stats/adegenet_output.RData`      

Open Rscript `01_scripts/adegenet_sep_pops.R` to import the Rdata from previous step.      
Analyzes the following contrasts:        
* Global with one representative from BC `global`
* All BC naturalized samples `bc_all`    
* BC spatial vs. temporal in BC (Hisnit and Serpentine) `bc_size`     
* All China `ch_all`
* BC, France, and China transfer to farm from nature `bc_sel`, `fr_sel`, `chr_sel`
* Domestication signatures for DPB vs PEN, ROS vs PIP, QDC vs CHN, RSC vs CHN, and GUR vs FRA    

Outputs with analysis-specific label to `11_adegenet_analysis`.     


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


## 6. Plotting outliers along the chromosome-level assembly
### A. Align loci against the reference genome
The populations run above (single snp) will output a fasta file with a consensus sequence per locus, entitled `09-diversity_stats/populations.loci.fa`. Align this against the most complete reference genome:      

Make a folder to work within:    
`mkdir 12_genome_plot`      

Move marker file to the folder:    
`cp 09-diversity_stats/populations.loci.fa 12_genome_plot`      

Do a bit of an update on the fasta file to make sure that all of your records are unique:    
`cat 12_genome_plot/populations.loci.fa | awk '{ print $1"-"$2 }' - | sed 's/\[//g' | sed 's/\,//g' | sed 's/-$//g' | grep -vE '^#' - > 12_genome_plot/populations.loci_fullname.fa`     

First download the assembly from http://dx.doi.org/10.12770/dbf64e8d-45dd-437f-b734-00b77606430a 
Citation: Gagnaire _et al._, 2018. Analysis of Genome-Wide Differentiation between Native and Introduced Populations of the Cupped Oysters _Crassostrea gigas_ and _Crassostrea angulata_, Genome Biology and Evolution, *10*(9), pp. 2518–2534, https://doi.org/10.1093/gbe/evy194 .     

Index the genome:    
`bwa index C_gigas_Assembly_1.fa.gz`      

Align all loci (single consensus locus per marker) against the reference genome:        
```
GENOME=/home/ben/Documents/genomes/C_gigas_Assembly_1.fa.gz ; bwa mem -t 6 $GENOME 12_genome_plot/populations.loci_fullname.fa > 12_genome_plot/populations.loci.sam
samtools view -Sb -F 4 -q 1 12_genome_plot/populations.loci.sam > 12_genome_plot/populations.loci.bam
```

Obtain alignment info per marker, saving the top mapq:      
`samtools view 12_genome_plot/populations.loci.bam | awk '{ print $1 "," $3 "," $4 }' - > 12_genome_plot/aligned_marker_info.csv`     


## Other Utilities 
Jump to:    
[Evaluate positions of SNPs in tags](#evaluate-positions-of-snps-in-tags)    
[Evaluate number of reads used in output](#evaluate-number-of-reads-used-in-output)    
[Run fastStructure](#run-faststructure)


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
