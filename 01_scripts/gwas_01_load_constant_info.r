# Plot differentiation statistics across genome
# Load constant info files
# inputs include:
# file "13_selection/chrom_pos_id.csv"
# output for "outlier_detection_rda.R" and 
# output for "outlier_detection_pcadapt.R" and 

#### Front Matter ####
# Clean space
# rm(list=ls())

# Load libraries
require("tidyr")
#install.packages("seqinr")
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

## Info
# sessionInfo()

### 1. Read in info files (Only the general info) ####
## Read in alignments (mname, chr, pos)
alignments <- read.csv(file = "12_genome_plot/aligned_marker_info.csv", header = F)
colnames(alignments) <- c("mname", "target.chr", "target.loc")
head(alignments)

# Alignment checking
length(unique(alignments$mname)) # Is there only a single alignment per input marker?
length(alignments$mname) # seems to be a handful of multiple alignments
# TODO# Keep top mapping locus (using mapq, keeping only one row per alignment query)

## Read in map file (interp with FST data)
map <- read.table(file = "11-adegenet_analysis/populations.plink.map")
map <- map[,c(1,2)] # Only need first two cols (stacks chr, locus_pos)
colnames(map) <- c("stacks.chr", "locus_pos")
head(map)
map$locus_pos <- as.character(map$locus_pos) # make mname character
map <- map[with(map, order(locus_pos)), ] # Order by locus_pos for data checking later
head(map)

## Read in genome file (for chr names and lengths)
# Provide full path to genome
genome.path <- "~/Documents/genomes/C_gigas_Assembly_1.fa"
rfsq <- read.fasta(file = genome.path, seqtype = "DNA") # read in genome
chr.names <- names(rfsq) # characterize chr names
chr.lengths <- as.numeric(getLength(object = rfsq)) # characterize chr length
chr.info <- as.data.frame(cbind(chr.names, chr.lengths)) # make object to contain this info
chr.info$chr.lengths <- as.numeric(as.character(chr.info$chr.lengths))
chr.info$chr.names <- as.character(chr.info$chr.names)

# Find cumulative position of chromosomes
chr.info$cum.lengths <- cumsum(chr.info$chr.lengths)
str(chr.info)

# Output is:
head(alignments)
head(map)
head(chr.info)


#### 1.1 Find cumulative pos of alignments ####
head(alignments)

alignments_w_cum_pos <- alignments
head(alignments_w_cum_pos)
colnames(alignments_w_cum_pos) <- c("mname","chr","pos") # rename columns

# Sort by chromosome, then position
alignments_w_cum_pos <- alignments_w_cum_pos[order(alignments_w_cum_pos$chr, alignments_w_cum_pos$pos), ]
head(alignments_w_cum_pos)

### subset rows that did not map
#unplaced.alignments_w_cum_pos <- alignments_w_cum_pos[which(is.na(alignments_w_cum_pos$chr)) , ]
#dim(unplaced.alignments_w_cum_pos) 
#dim(alignments_w_cum_pos)

# Keep only rows that mapped
#alignments_w_cum_pos <- alignments_w_cum_pos[which(!is.na(alignments_w_cum_pos$chr)) , ]
#dim(alignments_w_cum_pos)
#head(alignments_w_cum_pos)
#tail(alignments_w_cum_pos)

# Add a new column to alignments_w_cum_pos that is a continuous value of total position on the genome
cum.pos <- NULL; cum.pos.string <- NULL; chr <- NULL; identical.row <- NULL; adding.factor <- NULL

for(j in 1:nrow(alignments_w_cum_pos)){
  
  # Find the chromosome for this marker
  chr <- as.character(alignments_w_cum_pos$chr[j])
  print(chr)
  
  # Determine details of this chromosome in the chr.info file
  identical.row <- which(chr.info$chr.names==chr)
  
  # Identify the cumulative position of the preceeding chromosome (unless at first chr)
  if(identical.row > 1){
    adding.factor <- chr.info[identical.row - 1, "cum.lengths"]
  } else {
    adding.factor <- 0
  }
  
  # Add this cumulative position to the current marker
  cum.pos[j] <- alignments_w_cum_pos$pos[j] + adding.factor
}

# Add the cumulative position vector to the alignments_w_cum_pos
alignments_w_cum_pos <- cbind(alignments_w_cum_pos, cum.pos)
head(alignments_w_cum_pos)
tail(alignments_w_cum_pos)


#### 1.2 Find chr midpoints for plotting ####
position.of.chr.label <- NULL

for(k in 1:nrow(chr.info)){
  position.of.chr.label[k] <- (chr.info$chr.lengths[k] / 2) + chr.info$cum.lengths[k - 1]
}

position.of.chr.label[is.na(position.of.chr.label)] <- (chr.info$chr.lengths[1] / 2) # replace the first one due to issue of -1 from 1 position
position.of.chr.label

alignments <- alignments_w_cum_pos

# 
# Output is:
head(alignments)
head(map)
head(chr.info)

