# To reduce samples for going into fineRADstructure
# Clean space
# rm(list=ls())
setwd("~/Documents/01_moore_oyster_project/stacks_workflow_all_data/08-fineRADstructure") ## The directory where the files are located

# Input haplotype tsv file
haplotypes <- read.delim2("batch_1.haplotypes.tsv.fineRADpainter.lociFilt.samples30%missFilt.txt")
dim(haplotypes)
colnames(haplotypes)
haplotypes[1:3,]

# Transpose
haplotypes.t <- t(x = haplotypes)
dim(haplotypes.t)
haplotypes.t[1:5, 1:5]

# Understand the data a bit more
# Here are the populations that you have in the data
samples <- rownames(haplotypes.t)[3:length(rownames(haplotypes.t))]

library(tidyr)
separate(data = as.data.frame(samples), sep = "_")
samples.df <- as.data.frame(samples)

pop <- samples.df %>% separate(samples, c("pop","id"), sep = "_")
pop <- pop$pop
pops.all <- unique(pop)
pops.all

# Remove unneeded samples
reduced.haplotypes.t <- haplotypes.t[-grep("PENF|PIP|SER|HIS", rownames(haplotypes.t)), ]
dim(reduced.haplotypes.t)


# Transpose back
reduced.haplotypes <- t(reduced.haplotypes.t)

# Write out
write.table(x = reduced.haplotypes, file = "batch_1.haplotypes.tsv.fineRADpainter.lociFilt.samples30%missFilt.txt"
            , quote = F, sep = "\t"
            , row.names = F)
