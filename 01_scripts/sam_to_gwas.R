# Identify cumulative positions of markers on fasta

#### Front Matter ####
# Clean space
# rm(list=ls())

# Install libraries
source("http://bioconductor.org/biocLite.R")
# biocLite("ShortRead")
# biocLite("SimRAD")
# require("SimRAD")

#install.packages("seqinr")
require("seqinr")

# Set working directory
# Xavier
setwd("~/Documents/01_moore_oyster_project/stacks_workflow")

# Provide full path to genome
genome.path <- "~/Documents/genomes/C_gigas_assembly_1_all_chr.fa"

# # Read in genome
# rfsq <- ref.DNAseq(FASTA.file = genome.path, subselect.contigs = F)
# length(rfsq)

# Read in genome
rfsq <- read.fasta(file = genome.path, seqtype = "DNA")
str(rfsq)

# Find chr names and lengths
chr.names <- names(rfsq)
chr.lengths <- as.numeric(getLength(object = rfsq))

chr.info <- as.data.frame(cbind(chr.names, chr.lengths))
chr.info$chr.lengths <- as.numeric(as.character(chr.info$chr.lengths))
chr.info$chr.names <- as.character(chr.info$chr.names)
str(chr.info)

# find cumulative position of chromosomes
chr.info$cum.lengths <- cumsum(chr.info$chr.lengths)
str(chr.info)



#### Read in alignment info ####
marker.info <- read.csv(file = "11-other_stats/aligned_marker_info.csv", header = F, col.names = c("mname","chr","pos"))
marker.info

# remove unaligned
marker.info <- marker.info[marker.info$chr!="*", ]
dim(marker.info)

head(marker.info)

### WORKING HERE ###
# translate into continuous
cum.pos <- NULL
cum.pos.string <- NULL
chr <- NULL
identical.row <- NULL
adding.factor <- NULL


for(i in 1:nrow(marker.info)){
  chr <- marker.info$chr[i]
  chr <- as.character(chr)
  print(chr)
  identical.row <- which(chr.info$chr.names==chr)
  adding.factor <- chr.info[identical.row - 1, "cum.lengths"]
  cum.pos <- marker.info$pos + adding.factor
  cum.pos.string <- c(cum.pos.string, cum.pos)
  
  #cum.pos[i] <- marker.info$pos[i] + 
}


### WORKING HERE ###
