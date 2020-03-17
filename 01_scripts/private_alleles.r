# Separate Analysis of individual populations
# Requires the running of 'adegenet.R'

#### Front Matter ####
# Clean space
# rm(list=ls())

# Load required libraries
library("purrr")
library("dplyr")
library("ape")
library("pegas")
library("seqinr")
library("ggplot2")
library("devtools")
library("hierfstat")
library("adegenet")
library("lattice")
library("gplots")
#install.packages("poppr")
library("poppr")


# Set working directory for stark or xavier
# if on a different system, prompt for working directory
if(Sys.info()["nodename"] == "stark"){ 
  print("On Stark, ready to go")
  setwd("/mnt/data/bsuther/01_moore_oyster_project/stacks_workflow/") # stark
} else if(Sys.info()["nodename"] == "Xavier"){
  print("On Xavier, ready to go")
  setwd("~/Documents/01_moore_oyster_project/stacks_workflow_all_data/") # Xavier
} else {
  print("You are on an unrecognized system, please set working directory manually")
}

## Info
# sessionInfo()

# Set variables
output.dir <- "11-adegenet_analysis/"

#### 01. Input data and prepare ####
## Load inputs
# Load part 1 results
load(file = paste0(output.dir, "adegenet_output.RData"))
my.data.gid

# Load colours file
my_cols.df <- read.csv(file = "../ms_oyster_popgen/00_archive/my_cols.csv", stringsAsFactors = F)
str(my_cols.df)

#### 02. Private alleles

# Tabulate alleles the occur in only one population. 
private_alleles(gid = my.data.gid)

# Do separately for different contrasts
## Prepare individual analyses of genind by separating all pops
sep.obj <- seppop(x = my.data.gid)
names(sep.obj)

# Create list of various combinations of repooled pops
# "global" is the comparison of one representative from each main pop
datatype.list <- list()
datatype.list[["global"]] <- repool(sep.obj$PEN
                                    , sep.obj$CHN, sep.obj$QDC, sep.obj$RSC
                                    , sep.obj$DPB, sep.obj$ROS, sep.obj$GUR
                                    , sep.obj$FRA
                                    , sep.obj$JPN
)

datatype.list[["all"]] <- my.data.gid

datatype.list[["bare.minimum"]] <- repool(sep.obj$PEN
                                          , sep.obj$CHN
                                          , sep.obj$DPB
                                          , sep.obj$GUR
                                          
)


## Select the dataset
## Main analysis:
global.gid <- datatype.list[["all"]]

## Test analyses
#global.gid <- datatype.list[["bare.minimum"]]
#global.gid <- datatype.list[["global"]]

## Create a df that defines the strata for each individual in the rows
strata.df <- as.data.frame(as.character(pop(global.gid)), stringsAsFactors = FALSE)
colnames(strata.df)[1] <- "indiv.pop"
table(strata.df)
str(strata.df)
# unique(strata.df$repunit)

# Add repunit
strata.df$repunit <- strata.df$indiv.pop

# Grouping similar pops into the same strata (not all will be present in all datasets, but that is OK)
strata.df$repunit <- gsub(x = strata.df$repunit, pattern = "CHNF", replacement = "CHN")
strata.df$repunit <- gsub(x = strata.df$repunit, pattern = "FRAF", replacement = "FRA")
strata.df$repunit <- gsub(x = strata.df$repunit, pattern = "PENF", replacement = "BC")
strata.df$repunit <- gsub(x = strata.df$repunit, pattern = "PIP", replacement = "BC")
strata.df$repunit <- gsub(x = strata.df$repunit, pattern = "PEN", replacement = "BC")
strata.df$repunit <- gsub(x = strata.df$repunit, pattern = "SER", replacement = "BC")
strata.df$repunit <- gsub(x = strata.df$repunit, pattern = "HIS", replacement = "BC")

unique(strata.df$repunit)

# Add strata object to gid
strata(global.gid) <- strata.df
table(strata(global.gid))

## Characterize sample sizes
table(strata(global.gid, ~indiv.pop))
table(strata(global.gid, ~repunit, combine = FALSE))
table(strata(global.gid, ~repunit/indiv.pop, combine = FALSE))

# Now try to use strata
# # Tabulate alleles the occur in only one population. 
# private_alleles(gid = my.data.gid)

# find private alleles per population
per_pop.privallele <- private_alleles(global.gid, alleles ~ indiv.pop)
per_repunit.privallele <- private_alleles(global.gid, alleles ~ repunit)
per_repunit.privallele

# Observe output object
dim(per_repunit.privallele)
per_repunit.privallele[,1:5]

# Quantify the number of observations of each private allele in each population
# e.g.,
table(per_repunit.privallele["BC",])

# for all repunits
for(i in 1:nrow(per_repunit.privallele)){
  
  print(paste0("The repunit ", rownames(per_repunit.privallele)[i], " has the following observed priv alleles"))
  print(table(per_repunit.privallele[i, ]))
  
}

# Number of observed alleles per locus
# private_alleles(Pinf, locus ~ Country, count.alleles = TRUE)
# private_alleles(global.gid, count.alleles = TRUE) # I don't know the difference here..

# Get raw number of private alleles per locus
pal <- private_alleles(global.gid, locus ~ repunit, count.alleles = FALSE)
table(pal) # this only gives a 0 or 1, does not count alleles, allows one to see exact how many private alleles exist
# per repunit
# This shows the number of rows that have a private allele
rowSums(pal)


# Plot 
global.priv <- private_alleles(global.gid, alleles ~ repunit, report = "data.frame")
ggplot(global.priv) + geom_tile(aes(x = population, y = allele, fill = count))

