# Load genome alignments and the genome
# input(s) include: (see 'Set variables' below)

#### Front Matter ####
# Clean space
# rm(list=ls())

# Load libraries
#install.packages("seqinr")
require("tidyr")
require("seqinr")

# Set working directory for stark or xavier
# if on a different system, prompt for working directory
if(Sys.info()["nodename"] == "stark"){ 
  print("On Stark, ready to go")
  setwd("/mnt/data/01_moore_oyster_project/stacks_workflow/") # stark
} else if(Sys.info()["nodename"] == "Xavier"){
  print("On Xavier, ready to go")
  setwd("/hdd/01_moore_oyster_project/stacks_workflow/") # Xavier
} else {
  print("You are on an unrecognized system, please set working directory manually")
}

# Set variables
alignment_info.FN <- "12_genome_plot/aligned_marker_info.csv"
genome_assembly.FN <- "~/Documents/genomes/C_gigas_Assembly_1.fa"

### 1. Read in info files (Only the general info) ####
## Read in alignments (mname, chr, pos)
alignments <- read.csv(file = alignment_info.FN, header = F)
colnames(alignments) <- c("mname","chr","pos")
head(alignments)

# # Alignment checking
# length(unique(alignments$mname)) # Is there only a single alignment per input marker?
# length(alignments$mname) # seems to be a handful of multiple alignments
# # TODO# Keep top mapping locus (using mapq, keeping only one row per alignment query)

## Read in genome file (for chr names and lengths)
rfsq <- read.fasta(file = genome_assembly.FN, seqtype = "DNA")
chr.names <- names(rfsq)
chr.lengths <- as.numeric(getLength(object = rfsq))
chr.info <- as.data.frame(cbind(chr.names, chr.lengths)) # make object to contain this info
chr.info$chr.names   <- as.character(chr.info$chr.names)
chr.info$chr.lengths <- as.numeric(as.character(chr.info$chr.lengths))
chr.info$cum.lengths <- cumsum(chr.info$chr.lengths) # Find cumulative position of chromosomes
str(chr.info)

# Outputs:
head(alignments)
head(chr.info)
