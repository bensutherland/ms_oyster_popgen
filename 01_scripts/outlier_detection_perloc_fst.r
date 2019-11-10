# Outlier detection pcadapt but with perloc FST

#### Front Matter ####
# Clean space
# rm(list=ls())

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

## Install packages
# install.packages("vcfR")
library("vcfR")
library("adegenet")
library("hierfstat")

# Read in a VCF (need to loop probably)
my_vcf <- read.vcfR(file = "13_selection/DPB_PEN.recode.vcf")

# Convert VCF to genind
my_genind <- vcfR2genind(x = my_vcf) # convert to genind
my_genind$tab[1:5,1:5] # view data

# Assign pop attribute
indiv_names <- indNames(my_genind)
indiv_pops <- gsub(pattern = "_.*", replacement = "", x = indiv_names)
pop(my_genind) <- indiv_pops
pop(my_genind)

# View marker names
head(colnames(my_genind$tab), n = 25)

# Convert genind to hierfstat
my_hf <- genind2hierfstat(dat = my_genind, pop = pop(my_genind))
dim(my_hf)
my_hf[1:5,1:5]
colnames(my_hf) <- gsub(pattern = "X", replacement = "CLocus_", colnames(my_hf))
colnames(my_hf) <- gsub(pattern = "\\..*", replacement = "", colnames(my_hf))
# Should now match qvals

# Calculate per-locus FST
my_basic_stats <- basic.stats(data = my_hf, diploid = TRUE, digits = 4)
my_perloc <- my_basic_stats$perloc
my_perloc$mname <- rownames(my_perloc)
head(my_perloc)


# Let's look at the qvals file
my_qvals <- read.csv(file = "13_selection/qvals_DPB_PEN.csv")
head(my_qvals) # note mname should match the marker names in the genind
my_qvals$mname <- gsub(pattern = ">", replacement = "", my_qvals$mname)
my_qvals$mname <- gsub(pattern = "-.*", replacement = "", my_qvals$mname)


# if the following are equal, then all CLocus numbers are unique and we can match the hf and the qvals quite easily
length(my_qvals$mname) 
length(unique(my_qvals$mname))


# MERGE FST AND QVAL
my_FST_and_qval <- merge(x = my_perloc, y = my_qvals, by = "mname")
head(my_FST_and_qval)
# Drop NAs
table(is.na(my_FST_and_qval$qval))
my_FST_and_qval <- my_FST_and_qval[!is.na(my_FST_and_qval$qval),]  # Keep non-NA
dim(my_FST_and_qval)

plot(my_FST_and_qval[my_FST_and_qval$qval < 0.01, "Fst"] ~ my_FST_and_qval[my_FST_and_qval$qval < 0.01, "qval"]
     , xlab = "qval", ylab = "Fst"
     )
