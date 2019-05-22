# Originally was created to input markers out of HWP from stacks output through PLINK format
# Since the HWE filter is going to be used for all markers, this script is no longer active

#### Front Matter ####
# Clean space
# rm(list=ls())

## Install and load packages
library("hierfstat")
library("adegenet")
library("purrr")
library("dplyr")
library("ape")
library("pegas")
library("seqinr")
library("ggplot2")
library("devtools")
library("diveRsity")

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
output.dir <- "11-adegenet_analysis/"
input.FN <- "24-selection/populations_single_snp.raw"

#### 1. Import data ####
my.sel.data.gl <- read.PLINK(file = input.FN)
my.sel.data.gl

#### 2. Convert genlight to genind ####
# Convert genlight to matrix
my.data.mat <- as.matrix(my.sel.data.gl)
my.data.mat[1:5,1:5]

# Translate the number of minor allele to genind format
my.data.mat[my.data.mat == 0] <- "1/1" #homozygote reference
my.data.mat[my.data.mat == 1] <- "1/2" #heterozygote
my.data.mat[my.data.mat == 2] <- "2/2" #homozygote alternate
my.data.mat[1:5,1:5]

# Convert matrix to genind
my.sel.data.gid <- df2genind(my.data.mat, sep = "/", ploidy = 2) # convert df to genind
my.sel.data.gid

# Transfer pop attributes
pop(my.sel.data.gid) <- pop(my.sel.data.gl) 

# Save out
save(my.sel.data.gid, file = "24-selection/my_sel_data.Rdata")

# Then go to 'adegenet_sep_pops.R' for full analysis