#### Front Matter ####

# Clean space
# rm(list=ls())

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

#### Import Genomics (genlight) ####
#getClassDef("genlight") # info on the class type
my.data <- read.PLINK(file = "06-stacks_rx/batch_1.raw"
                      #, map.file = "06-stacks_rx/batch_1.plink.map"
                      )
my.data

# Get information about the data
nInd(my.data) # how many individuals?
nLoc(my.data) # how many loci?
indNames(my.data) # see names of individuals
nPop(my.data) # how many pops?
pop(my.data) # what are the pop names?
head(locNames(my.data, withAlleles = T), n = 20) # see allele names


# Convert data to a list if necessary
str(as.list(my.data))
str(as.matrix(my.data))
# and these can be subset this way as well

# Positional info
# temp <- density(x = position(my.data), bw=10) # there is no position info, for some reason
# my.data


# Basic analyses
# genlight plot to visualize the number of second alleles across individuals and SNPs
glPlot(x = my.data, posi="topleft") 

# # Compute the sum of second alleles for each SNP
# glSum.my.data <- glSum(my.data)
# plot(glSum.my.data[1:1000])
# 
# # compute the number of missing values per locus
# glNA.my.data <- glNA(my.data)
# plot(glNA.my.data[1:1000])
# # and others, see tutorial (glMean, glVar, glDotProd)

# Create density plot of minor allele frequencies
myFreq <- glMean(my.data)
hist(myFreq, proba=T, col="gold", xlab = "Allele frequencies"
     , main = "Distribution of (second) allele frequencies"
     , ylim = c(0,8)
     )
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

#### Build neighbour-joining tree ####
# Determine distance between individuals
# First use seploc to create ten blocks of loci, each to be computed in parallel
x <- seploc(my.data, n.block=10, parallel=F)
names(x); x[1:2] # can query the blocks

## Use dist in a lapply loop to compute pairwise distances between individuals for each block
lD <- lapply(x, function(e) dist(as.matrix(e)))
class(lD)
names(lD)
class(lD[[1]]) # lD is a list of distance matrices bw pairs of indiv

# Obtain general distance matrix by summing these
D <- Reduce("+", lD)

# Now carry on to make a neighbor-joining tree using ape
plot(nj(D), type="fan")


#### Run a Principle Components Analysis on the SNP data ####
pca1 <- glPca(my.data)
scatter(pca1, posi = "topleft")
title("PCA of the 4 PCS")

scatter(x = pca1, posi = "topleft", xax = 1, yax = 2) # show PC1 and PC2
scatter(x = pca1, posi = "topleft", xax = 3, yax = 4) # show PC3 and PC4

# Plot the loading values of the different markers into the PCA
loadingplot(x = pca1)

# Plot a colorplot using the PC scores of the samples
#basically colors samples in Red Green Blue where colors correspond to loadings of PC1, PC2, PC3
myCol <- colorplot(pca1$scores, pca1$scores, transp=T, cex=4)
abline(h=0,v=0, col = "grey")
# add.scatter.eig(pca1$eig[1:40], 2,1,2, posi="topleft", inset=0.05, ratio=0.3) # adds eigenplot to image
#The colors indicate the values for PC1 (Red), PC2 (Green), PC3 (Blue)
#The nice thing is that it brings in a third dimension..

# Can bring these colors into the Neighbour-joining tree
plot(nj(D), type="fan", show.tip=F)
tiplabels(pch=20, col=myCol, cex=4)
title("NJ tree of my.data")


#### Discriminant Analysis of Principal Components (DAPC) ###
# DAPC implemented for genlight by an approp. method for find.clusters and dapc generics
# Can run find.clusters and dapc on genlight objects and the approp. functions will be used
# 1st transforms data (possibly), submits to PCA, 2nd PCs are submitted to Linear Discriminant Analysis (LDA)

dapc1 <- dapc(my.data, n.pca = 10, n.da = 1) 
# n.pca = number axes to be retained; n.da = number of axes retained in Discriminant Analysis step
dapc1 # variance = 0.17

# Density plot of samples along discriminant function 1
scatter(dapc1, scree.da = F, bg = "white", posi.pca = "topright", legend = T
        , txt.leg=c("Pipestem","Pendrell","DeepBay","Rosewall","BayneSpat","PendrellFarm")
        )

# Composition plot
compoplot(dapc1
          #, lab="" # comment out if you want sample labels
          , txt.leg=c("Pipestem","Pendrell","DeepBay","Rosewall","BayneSpat","PendrellFarm")
            )

# Loading plot
loadingplot(dapc1$var.contr, thres=1e-3)



