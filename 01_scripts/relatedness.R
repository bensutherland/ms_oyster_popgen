# related for pairwise comparison

#### 00. Front Matter ####
# Clean space
# rm(list=ls())

## Install packages:
# install.packages("related", repos="http://R-Forge.R-project.org")
# install.packages("Demerelate")
# install.packages("BiocManager")
# require("BiocManager")
# BiocManager::install("SNPRelate")
# BiocManager::install("qvalue")
# install.packages("dartR") # To prepare a genlight file to go into demerelate
# install.packages("tidyr")
 
require("Demerelate")
require("related")
require("SNPRelate")
require("qvalue")
require("dartR")
require("tidyr")

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
input.dir <- "11-adegenet_analysis/"
output.dir <- "09-diversity_stats/"

#### 01. Input data and prepare
# Load part 1 results
load(file = paste0(input.dir, "adegenet_output.RData"))

# Your genlight data:
my.data

# Create demerelate df from genlight object
my.data.demerelate <- gl2demerelate(gl = my.data)
my.data.demerelate[1:5,1:5]


#### 02. Relatedness with demerelate ####
#### Does this work?  ####
# dem_results <- Demerelate(inputdata = my.data.demerelate[,1:60], value = "Mxy"
#            , file.output = F, object = T, pairs = 10
#            #, 
#            )


#### 03. Relatedness with related ####
# Info on related input can be obtained here: https://rdrr.io/rforge/related/man/related-package.html
## genotype file should have one column of individual identifiers and 
## then 2 columns for each locus (one for each allele). No other columns are allowed. 
## Missing data should be formatted as zeros ("0").
## The file should NOT contain a header row.

# Convert to related input
str(my.data.demerelate[1:5,1:5])

# Remove pop column
my.data.related <- my.data.demerelate[,-2] # drop pop column
str(my.data.related[1:10,1:10])

# Convert NA to "0" for related format
my.data.related[is.na(my.data.related)] <- 0 
str(my.data.related[1:10,1:10])
my.data.related[1:5,1:5]

# Save out as text file to be able to read in via readgenotypedata
write.table(x = my.data.related, file = "11-adegenet_analysis/my_data_related.txt"
            , quote = F, sep = "\t", row.names = F, col.names = F)

# Read in via readgenotypedata
my.data.related <- readgenotypedata(genotype.data = "11-adegenet_analysis/my_data_related.txt")

names(my.data.related)

# coancestry analysis w/ related
output <- coancestry(genotype.data = my.data.related$gdata
                     , lynchrd = 2
                     , quellergt = 2
                     , wang = 2
                     ) # all settings default from website

# Save out results
date <- format(Sys.time(), "%Y-%m-%d")
save.image(file = paste0(output.dir, "kinship_analysis_", date, ".Rdata"))
# Load results back up (#to build)




## Values: 
# relatedness: df w/ all pairwise est. of relatedness
# 1. integer for pair #; 2. indiv 1 ID; 3. indiv 2 ID; 4. group assignment; 5-11. relatedness estim.
# delta7 & delta8: df w/ delta7 & delta8 est. for relatedness estim. that use it
# inbreeding: df w/ inbreeding est. per indiv., as used in relatedness est (one row per indiv)

# See what variables are available for relatedness output
head(output$relatedness, n = 5)
colnames(output$relatedness)

# Identify which samples will be same-on-same
comp.names <- unique(output$relatedness$group)
comp.names
comp.names <- as.data.frame(comp.names) # make df for sep below
comp.names
comp.names <- as.data.frame(sort(comp.names$comp.names)) # sort as boxplot will
colnames(comp.names) <- "comp.names"
head(comp.names)

test <- separate(data = comp.names, col = comp.names, into = c("first", "second"), sep = 2) # split
head(test)

test$TF <- test$first==test$second
head(test)

for(i in 1:nrow(test)){
  if(test$TF[i]==TRUE){
    test$colour[i] <- "red"
  }else {test$colour[i] <- "white"}
}

head(test)

colours <- as.character(test$colour)


## Set up a wrapper to plot all types
datatypes <- c("wang", "ritland","quellergt")
# datatype <- "wang"
# datatype <- "ritland"
# datatype <- "quellergt"

## Loop
for(i in 1:length(datatypes)){
  datatype <- datatypes[i]
  
  pdf(file = paste0("09-diversity_stats/", "relatedness_", datatype, "_", date, ".pdf"), width = 12, height = 8)
  boxplot(output$relatedness[[datatype]] ~ output$relatedness$group
          , col = colours
          , las = 2
  )
  abline(h = 0, lty = 2)
  dev.off()
  
}
