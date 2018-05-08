#### Front Matter ####

# Clean space
# rm(list=ls())

#par(mfrow=c())

# Set working directory
setwd("~/Documents/01_moore_oyster_project/stacks_workflow_2018-05-03/")

# Resource packages required:
#install.packages("purrr")
library("purrr")
#install.packages("dplyr")
library("dplyr")

# install.packages("hierfstat")
library(hierfstat)

# Other useful packages
library("ape")
library("pegas")
library("seqinr")
library("ggplot2")

# Install packages
library("adegenet")

# # Vignettilable 
# vignette("hierfstat")

#### import ####

# Import a genpop file (output from stacks) into a genind file  
my.genind <- import2genind(file = "06-stacks_rx/batch_1.gen")

## Evaluate genpop file ##
nInd(my.genind)
nLoc(my.genind)

indNames(my.genind)

nPop(my.genind)
pop(my.genind) # apparently this is not properly set yet
pop(my.genind) <- c(rep("Pendrell", times = 25), rep("Pipestem", times = 23))


head(locNames(my.genind, withAlleles = T) , n = 20) # retrieves loci and allele names



#### Genomics ####
getClassDef("genlight")

my.data <- read.PLINK(file = "07-filtered_vcfs/plink.raw"
                      #, map.file = "06-stacks_rx/batch_1.plink.map"
                      )

my.data

# this can be converted into a list or matrix as follows
str(as.list(my.data))
str(as.matrix(my.data))
# and these can be subset this way as well.

# Positional info
temp <- density(x = position(my.data), bw=10) # there is no position info, for some reason
my.data


# Basic analyses
# genlight plot to visualize the number of second alleles across individuals and SNPs
glPlot(x = my.data, posi="topleft") 

# Compute the sum of second alleles for each SNP
glSum.my.data <- glSum(my.data)
plot(glSum.my.data[1:1000])

# compute the number of missing values per locus
glNA.my.data <- glNA(my.data)
plot(glNA.my.data[1:1000])

# and others, see tutorial (glMean, glVar, glDotProd)
# Create density plot of minor allele frequencies
myFreq <- glMean(my.data)
hist(myFreq, proba=T, col="gold", xlab = "Allele frequencies"
     , main = "Distribution of (second) allele frequencies")
temp <- density(myFreq)
lines(temp$x, temp$y, lwd=3)


# not working?
# # Positions of NAs in the alignments
# head(glNA(my.data), 20)
# temp <- density(glNA(my.data), bw =10)
# plot(temp, type="n", xlab = "Position in the alignment"
#      , main = "Location of the missing value"
#      )
# polygon(c(temp$x, rev(temp$x)), c(temp$y, rep(0,length(temp$x)))
#         , col=transp("blue", 0.3))
# points(glNA(my.data), rep(0, nLoc(my.data)), pch="|", col="blue")
# #

### Determine distance between individuals
# use seploc to create ten blocks, each to be computed in parallel
x <- seploc(my.data, n.block=10, parallel=F)
names(x)
x[1:2]

# use dist in a lapply loop to compute pairwise distances between individuals for each block
lD <- lapply(x, function(e) dist(as.matrix(e)))
class(lD)
names(lD)
class(lD[[1]])

# lD is a list of distance matrices b/w pairs of indivs. 
# obtain general distance matrix by summing these
D <- Reduce("+", lD)

# Now carry on to make a neighbor-joining tree using ape
library(ape)

plot(nj(D), type="fan")








pca1 <- glPca(my.data)
scatter(pca1)













