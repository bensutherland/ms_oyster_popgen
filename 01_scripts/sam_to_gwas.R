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
setwd("~/Documents/01_moore_oyster_project/stacks_workflow_no_Qc/")

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
dim(marker.info)
head(marker.info)

# Sort by chromosome, then position
marker.info <- marker.info[order(marker.info$chr, marker.info$pos), ]
tail(marker.info)

# TODO: should also sort the chr.info file

# Add a new column to marker.info that is a continuous value of total position on the genome
cum.pos <- NULL; cum.pos.string <- NULL; chr <- NULL; identical.row <- NULL; adding.factor <- NULL

for(i in 1:nrow(marker.info)){
  # Find the chromosome for this marker
  chr <- as.character(marker.info$chr[i])
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
  cum.pos[i] <- marker.info$pos[i] + adding.factor
}


marker.info <- cbind(marker.info, cum.pos)
head(marker.info)
tail(marker.info)

# Find the midpoint b/w chromosomes for plotting purposes
position.of.chr.label <- NULL
for(i in 1:nrow(chr.info)){
  position.of.chr.label[i] <- (chr.info$chr.lengths[i] / 2) + chr.info$cum.lengths[i - 1]
}
position.of.chr.label[is.na(position.of.chr.label)] <- (chr.info$chr.lengths[1] / 2) # replace the first one due to issue of -1 from 1 position
position.of.chr.label



#### BRING IN DATA ####
datatypes <- c("ch.gid", "fr.gid", "bc.sel.gid")

pdf(file = "11-other_stats/selection_fst.pdf", width = 7, height = 12)
par(mfrow=c(3,1))

# loop to try all three
for(i in 1:length(datatypes)){
  datatype <- datatypes[i]

filename <- paste("11-other_stats/", datatype, "_perloc_stats.csv", sep = "")

# Matching names, marker info
match1 <- gsub(pattern = "CLocus_", replacement = "", x = marker.info$mname)
match2 <- gsub(pattern = "\\_.*", "", x = match1)
marker.info$mname <- match2
head(marker.info)

# Selected data
perloc_stats <- read.csv(file = filename)
head(perloc_stats)

# matching names, perloc stats
match3 <- gsub(pattern = "X", replacement = "", x = perloc_stats$X)
match4 <- gsub(pattern = "\\_.*", "", x = match3)
perloc_stats$X <- match4 
colnames(perloc_stats)[1] <- "mname"

head(perloc_stats)

# Merge together the data
positional.and.fst.data <- merge(x = marker.info, y = perloc_stats, by = "mname")
positional.and.fst.data <- positional.and.fst.data[order(positional.and.fst.data$chr, positional.and.fst.data$pos), ]
head(positional.and.fst.data)
dim(positional.and.fst.data)

#Plot
plot(x = positional.and.fst.data$cum.pos, y = positional.and.fst.data$Fst
    , xlab = "Cumulative Position (bp)", ylab = "Fst"
    , las = 1, ylim = c(-0.2, 0.4))
abline(v = chr.info$cum.lengths)
text(x = position.of.chr.label, y = 0.4, labels = chr.info$chr.names)
text(x = chr.info[nrow(chr.info), "cum.lengths"] / 2 , y = 0.35, labels = datatype)
abline(h = 0.15, lty = 2)
}

dev.off()



#### Directly compare data without any genome information ####
datatypes <- c("ch.gid", "fr.gid", "bc.sel.gid")
perloc_stats_list <- list()

# Loop
for(i in 1:length(datatypes)){
  datatype <- datatypes[i]
  print(datatype)
  filename <- paste("11-other_stats/", datatype, "_perloc_stats.csv", sep = "")

  # Selected data
  perloc_stats <- read.csv(file = filename)
  perloc_stats_list[[datatype]]  <- perloc_stats
}

str(perloc_stats_list)

pdf(file = "11-other_stats/selection_expt_all_by_all.pdf", width = 8, height = 8)
par(mfrow=c(2,2))

# France vs China
plot(perloc_stats_list[["ch.gid"]]["Fst"]$Fst ~ perloc_stats_list[["fr.gid"]]["Fst"]$Fst
     , ylim = c(0, 0.3)
     , xlim = c(0, 0.3)
     , ylab = "China (Fst)", xlab = "France (Fst)")
abline(h = 0.1, lty = 3)
abline(v = 0.1, lty = 3)

# BC vs France
plot(perloc_stats_list[["bc.sel.gid"]]["Fst"]$Fst ~ perloc_stats_list[["fr.gid"]]["Fst"]$Fst
     , ylim = c(0, 0.3)
     , xlim = c(0, 0.3)
     , ylab = "BC (Fst)", xlab = "France (Fst)")
abline(h = 0.1, lty = 3)
abline(v = 0.1, lty = 3)

# BC vs China
plot(perloc_stats_list[["bc.sel.gid"]]["Fst"]$Fst ~ perloc_stats_list[["ch.gid"]]["Fst"]$Fst
     , ylim = c(0, 0.3)
     , xlim = c(0, 0.3)
     , ylab = "BC (Fst)", xlab = "China (Fst)")
abline(h = 0.1, lty = 3)
abline(v = 0.1, lty = 3)
dev.off()

# Find out which markers correspond to this:
perloc_ch.gid <- as.data.frame(perloc_stats_list$ch.gid)
perloc_fr.gid <- as.data.frame(perloc_stats_list$fr.gid)
perloc_bc.sel.gid <- as.data.frame(perloc_stats_list$bc.sel.gid)

FR <- perloc_fr.gid[which(perloc_fr.gid$Fst > 0.1), "X"]
BC <- perloc_bc.sel.gid[which(perloc_bc.sel.gid$Fst > 0.1), "X"]

intersect(x = FR, y = BC) # tells you which markers




