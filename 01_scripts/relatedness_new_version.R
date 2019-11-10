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
  setwd("/mnt/data/01_moore_oyster_project/stacks_workflow/") # stark
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


#### 02. Relatedness with related ####
# Info on related input can be obtained here: https://rdrr.io/rforge/related/man/related-package.html
## genotype file should have one column of individual identifiers and 
## then 2 columns for each locus (one for each allele). No other columns are allowed. 
## Missing data should be formatted as zeros ("0").
## The file should NOT contain a header row.

# Convert to related input
my.data.demerelate[1:5,1:5]

# Remove population column
my.data.related <- my.data.demerelate[,-2] # drop pop column
my.data.related[1:5,1:5]

# Convert NA to "0" for related format
my.data.related[is.na(my.data.related)] <- 0 
my.data.related[1:5,1:5]

# Write out as text file to be able to read in via readgenotypedata
write.table(x = my.data.related, file = "11-adegenet_analysis/my_data_related.txt"
            , quote = F, sep = "\t", row.names = F, col.names = F)

# Read in via readgenotypedata
my.data.related <- readgenotypedata(genotype.data = "11-adegenet_analysis/my_data_related.txt")

names(my.data.related)

### Need to change the names
my.data.related$gdata[1:5,1:5]

pop_name.df <- my.data.related$gdata$V1
head(pop_name.df)
str(pop_name.df)
pop_name.df <- as.data.frame(pop_name.df, stringsAsFactors = FALSE)
str(pop_name.df)
pop_name.df$pop_name <- as.character(pop_name.df$pop_name)
head(pop_name)

# Split the pop
temp <- separate(data = pop_name.df, col = pop_name, sep = "_", into = c("pop","ind"))

# Merge with stock code file
codes <- read.delim(file = "../ms_oyster_popgen/00_archive/pop_codes.txt", header = F, sep = "\t")
colnames(codes) <- c("pop", "pop_code")
codes$pop <- as.character(codes$pop)
str(codes)

head(temp)

# Re-attach this instead of the current ind name
temp_all.df <- merge(x = temp, y = codes, by = "pop", all.x = T)
head(temp_all.df)

# make a new vector with the pop code and ind
temp_all.df$new_ID <- paste0(temp_all.df$pop_code, "_", temp_all.df$ind)
head(temp_all.df)

# Replace original ID with this new one
my.data.related$gdata$V1 <- temp_all.df$new_ID

# Note that this is now using the stock code instead of the stock name
# , as the first two digits only will be retained
my.data.related$gdata[1:5,1:5]



# coancestry analysis w/ related
output <- coancestry(genotype.data = my.data.related$gdata
                     , lynchrd = 2
                     , quellergt = 2
                     , wang = 2
                     ) # all settings default from website

# Save out results
date <- format(Sys.time(), "%Y-%m-%d")
save.image(file = paste0(output.dir, "kinship_analysis_", date, ".Rdata"))


#### MANUAL LOAD BACK IN DATA #####
# Load results back up manually
load(file = "11-adegenet_analysis/kinship_analysis_2019-08-25.Rdata" )

## Values: 
# relatedness: df w/ all pairwise est. of relatedness
# 1. integer for pair #; 2. indiv 1 ID; 3. indiv 2 ID; 4. group assignment; 5-11. relatedness estim.
# delta7 & delta8: df w/ delta7 & delta8 est. for relatedness estim. that use it
# inbreeding: df w/ inbreeding est. per indiv., as used in relatedness est (one row per indiv)

# See what variables are available for relatedness output
head(output$relatedness, n = 5)
colnames(output$relatedness)

####
### Here forward: TAKEN FROM SIMPLE POP STATS ####
####
## The results are in the following formats
# relatedness: df w/ all pairwise est. of relatedness
# 1. integer for pair #; 2. indiv 1 ID; 3. indiv 2 ID; 4. group assignment; 5-11. relatedness estim.
# delta7 & delta8: df w/ delta7 & delta8 est. for relatedness estim. that use it
# inbreeding: df w/ inbreeding est. per indiv., as used in relatedness est (one row per indiv)

# See what variables are available for relatedness output
head(output$relatedness, n = 5)

# Extract the necessary parts
print("Extracting relevant sections of output data")
rel.wang <- output$relatedness[["wang"]]
rel.ritland <- output$relatedness[["ritland"]]
rel.quellergt <- output$relatedness[["quellergt"]]

group <- output$relatedness$group

# Make into df
output.df <- as.data.frame(
  cbind(  rel.wang
          , rel.ritland
          , rel.quellergt
          , group
  )
  , stringsAsFactors = F)

# Format sections of df
output.df$rel.wang <- as.numeric(output.df$rel.wang)
output.df$rel.ritland <- as.numeric(output.df$rel.ritland)
output.df$rel.quellergt <- as.numeric(output.df$rel.quellergt)

str(output.df)

# Make two vectors showing the pops in the contrast
output.df <- separate(data = output.df, col = group, into = c("pop1", "pop2"), sep = 2, remove = F) # split
str(output.df)

# compare pop1 and pop2 to identify if same
output.df$same <- output.df$pop1==output.df$pop2

head(output.df) # e.g. true
output.df[100:110,] # e.g. false
dim(output.df)

# Set variables
same_pops <- TRUE
datatype <- "SNP"


# Are all pop comparisons to be kept, or only same-on-same?
if(same_pops==TRUE){
  
  # Only keep same-on-same comparisons
  print("Keeping only same-on-same comparisons")
  output.df <- output.df[output.df$same==TRUE,]
  
}else if(same_pops==FALSE){
  print("Keeping all population comparisons")
}

dim(output.df)

## Set up a wrapper to plot all types
relatedness_metrics <- c("wang", "ritland","quellergt")

# Report
print("Plotting relatedness, output will be in 03_results")

#### DIFFERENT HERE ###

# Below, had to change to output.dir
#### END DIFFERENT HERE ###


## Loop to plot
for(i in 1:length(relatedness_metrics)){
  metric <- relatedness_metrics[i]
  
  pdf(file = paste0(output.dir, "relatedness_", metric, "_", date, ".pdf"), width = 20, height = 8)
  par(mar=c(7,3,3,3))
  boxplot(output.df[, paste0("rel.", metric)] ~ output.df$group
          #, col = comp_names.df$colour
          , las = 2
          , ylab = paste0("relatedness (", metric, ")")
  )
  abline(h = 0, lty = 2)
  dev.off()
  
}



#### ORIGINAL CODE #####

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
  
  pdf(file = paste0(output.dir, "relatedness_", datatype, "_", date, ".pdf"), width = 12, height = 8)
  boxplot(output$relatedness[[datatype]] ~ output$relatedness$group
          , col = colours
          , las = 2
  )
  abline(h = 0, lty = 2)
  dev.off()
  
}
