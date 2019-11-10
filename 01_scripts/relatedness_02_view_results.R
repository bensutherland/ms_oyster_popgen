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

# Make a new vector to not use a 2-letter grouping name, but rather the ind name
output$relatedness$group<- paste0(
        # Remove individual ID, keep pop ID
        gsub(x = output$relatedness$ind1.id, pattern = "_.*", replacement = "") 
        , "_"
        , gsub(x = output$relatedness$ind2.id, pattern = "_.*", replacement = "")
)

head(output$relatedness, n = 5)

# Extract the necessary parts
print("Extracting relevant sections of output data")
rel.wang <- output$relatedness[["wang"]]
rel.ritland <- output$relatedness[["ritland"]]
rel.quellergt <- output$relatedness[["quellergt"]]

group <- output$relatedness[["group"]]

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

#### Keep only same-on-same
# Make two vectors showing the pops in the contrast
output.df <- separate(data = output.df, col = group, into = c("pop1", "pop2"), sep = "_", remove = F) # split
str(output.df)

# compare pop1 and pop2 to identify if same
output.df$same <- output.df$pop1==output.df$pop2

head(output.df) # e.g. true
output.df[100:110,] # e.g. false
dim(output.df)

# Only keep same-on-same comparisons
print("Keeping only same-on-same comparisons")
output.df <- output.df[output.df$same==TRUE,]
dim(output.df)

### Rename pops simpler 
output.df$group <- gsub(x = output.df$group, pattern = "_.*", replacement = "")

### BRING IN COLOURS
## Set colors for datatypes
my.cols <- read.csv(file = "./../ms_oyster_popgen/00_archive/my_cols.csv")
rownames(my.cols) <- my.cols$my.pops
cols <- as.character(my.cols[c("PEN", "PENF", "PIP", "HIS", "SER"
                               , "DPB", "ROS", "FRA", "FRAF", "GUR"
                               , "QDC", "RSC", "CHN", "CHNF" ), 2]) # select colors into the same order that they are plotted

# ## Set up a wrapper to plot all types
relatedness_metrics <- c("wang", "ritland","quellergt")


# Create order of df
output.df$group <- factor(output.df$group , levels=c("PEN", "PENF", "PIP", "HIS", "SER"
                                                     , "DPB", "ROS", "FRA", "FRAF", "GUR"
                                                     , "QDC", "RSC", "CHN", "CHNF" ))


# To set wang for now..
i <- 1

# ## Loop to plot
# for(i in 1:length(relatedness_metrics)){
  metric <- relatedness_metrics[i]
  
  
  # test stats
  mod1 <- aov(output.df[, paste0("rel.", metric)] ~ output.df$group)
  summary(mod1)
  TukeyHSD(mod1)
  
  pdf(file = paste0(output.dir, "relatedness_", metric, "_", date, ".pdf"), width = 7, height = 5)
  par(mar=c(5,4,3,3))
  boxplot(output.df[, paste0("rel.", metric)] ~ output.df$group
          #, col = comp_names.df$colour
          , las = 2
          , ylab = paste0("relatedness (", metric, ")")
          , ylim = c(-0.4, 0.9)
          , col = cols
  )
  abline(h = 0, lty = 2)
  
  sig.indicators <- c("a", "b", "ad", "cdf", "cdf", "e", "adf", "b", "b", "adf", "f", "a", "ac", "ac")
  text(x = seq(1:length(levels(output.df$group))), y = 0.8
       , labels = sig.indicators
  )
  
  
  dev.off()
  
  
# }
