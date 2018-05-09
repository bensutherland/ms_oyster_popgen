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
library("ape")
library("pegas")
library("seqinr")
library("ggplot2")
library("adegenet")

# # Vignettilable 
# vignette("hierfstat")

#### 1. Import Genomics (genlight) ####
#getClassDef("genlight") # info on the class type
my.data <- read.PLINK(file = "06-stacks_rx/batch_1.raw"
                      #, map.file = "06-stacks_rx/batch_1.plink.map"
                      )
my.data

## Explore data attributes
nInd(my.data) # how many individuals?
nLoc(my.data) # how many loci?
indNames(my.data) # see names of individuals
nPop(my.data) # how many pops?
pop(my.data) # what are the pop names?
head(locNames(my.data, withAlleles = T), n = 20) # see allele names

## If you need to, you can convert this to a list or a matrix
#str(as.list(my.data))
#str(as.matrix(my.data))

## not currently applied: Positional Info ##
# temp <- density(x = position(my.data), bw=10) # there is no position info, for some reason
# my.data

#### 2. Basic Analysis (adegenet) ####
# genlight plot to visualize the number of second alleles across individuals and SNPs
glPlot(x = my.data, posi="topleft") 

# Create density plot of minor allele frequencies
myFreq <- glMean(my.data)
hist(myFreq, proba=T, col="gold", xlab = "Allele frequencies"
     , main = "Distribution of (second) allele frequencies"
     , ylim = c(0,8)
     )
text(x = 0.4, y = 7, labels = paste(nLoc(my.data), " loci", sep = "" ))
temp <- density(myFreq)
lines(temp$x, temp$y, lwd=3)


#### 3. Neighbour-joining tree, all loci ####
# Determine distance between individuals
# First use seploc to create a list of ten blocks of loci, each to be computed in parallel
x <- seploc(my.data, n.block=10, parallel=F)

## Use dist in a lapply loop to compute pairwise distances between individuals for each block
lD <- lapply(x, function(e) dist(as.matrix(e)))
class(lD)
names(lD)
class(lD[[1]]) # lD is a list of distance matrices bw pairs of indiv

# Obtain general distance matrix by summing these
D <- Reduce("+", lD)

# Plot a neighbor-joining tree using ape's nj function on the distance matrix
plot(nj(D), type="fan")


#### 4. Principle Components Analysis ####
pca1 <- glPca(my.data)

# Plot PCA
par(mfrow=c(2,1))
scatter(x = pca1, posi = "topleft", xax = 1, yax = 2)
#title("PC1 and PC2")
scatter(x = pca1, posi = "bottomleft", xax = 2, yax = 3) # show PC3 and PC4
#title("PC2 and PC3")

# Plot allele loadings
par(mfrow=c(1,1))
# Plot the loading values of the different markers into the PCA
loadingplot(x = pca1)

# Plot colorplot using the PC scores of the samples
myCol <- colorplot(pca1$scores, pca1$scores, transp=T, cex=4)
abline(h=0,v=0, col = "grey")
#The colors indicate the values for PC1 (Red), PC2 (Green), PC3 (Blue)

# Bring colorplot colors into the samples of the Neighbour-joining tree
plot(nj(D), type="fan", show.tip=F)
tiplabels(pch=20, col=myCol, cex=4)
title("NJ tree of my.data")

#### 5. Discriminant Analysis of Principal Components (DAPC) ####
# DAPC implemented for genlight by an approp. method for find.clusters and dapc generics
# Can run find.clusters and dapc on genlight objects and the approp. functions will be used
# 1st transforms data (possibly), submits to PCA, 2nd PCs are submitted to Linear Discriminant Analysis (LDA)

dapc1 <- dapc(my.data, n.pca = 10, n.da = 1) # n.pca = number axes to be retained; n.da = number of axes retained in Discriminant Analysis step
dapc1 # variance = 0.16

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


####6. Convert genlight to genind object ####
my.data.mat <- as.matrix(my.data)
my.data.mat[1:5,1:5]
my.data.mat[my.data.mat == 0] <- "1/1" #homozygote reference
my.data.mat[my.data.mat == 1] <- "1/2" #heterozygote
my.data.mat[my.data.mat == 2] <- "2/2" #homozygote alternate
#my.data.mat[my.data.mat == "NA"] <- "NA" # NA
my.data.mat[1:5,1:5]
my.data.gid <- df2genind(my.data.mat, sep = "/", ploidy = 2) # convert df to genind

# Generate genind population IDs
indNames(my.data.gid)
pop(my.data.gid) # currently null
pop(my.data.gid) <- gsub(x = indNames(my.data.gid), pattern = "\\_.*", replacement = "") # create pop ID
pop(my.data.gid)

# How many alleles?
table(nAll(my.data.gid))

# Descriptive Statistics on genind file
div <- summary(my.data.gid) # genetic diversity w/ adegenet
names(div) # contains

par(mfrow=c(1,3))
plot(div$Hobs
     , xlab="Loci number", ylab="Obs. Het."
     , main = "Observed heterozygosity per locus"
     , ylim = c(0,1))
# Note that stacks calculates obs het diff than other programs.. (missing data is not considered properly)
plot(div$Hexp
     , xlab="Loci number", ylab="Exp. Het."
     , main = "Expected heterozygosity per locus"
     , ylim = c(0,1))

# Plot as a function of each other:
plot(div$Hobs, div$Hexp
     , xlab = "Obs. Het.", ylab="Exp. Het."
     , main = "Expected Heterozygosity vs. Observed"
     , ylim = c(0,1))

# test : H0: Hexp = Hobs (using Bartlett test of homogeneity of variances)
exp.v.obs.test <- bartlett.test(list(div$Hexp, div$Hobs))
exp.v.obs.test
text(x = 0.6, y = 0.9, labels = paste("pval =", exp.v.obs.test$p.value, sep = ""))

##### 7. enter HIERFSTAT ####
# Create a hierfstat object (nloci+1 = columns; ninds = rows)
my.data.hfs <- genind2hierfstat(my.data.gid)

# Name the rows with the sample IDs
rownames(my.data.hfs) <- indNames(my.data.gid)

# OPTIONAL TO LIMIT DATA
#my.data.mini <- my.data.hfs[, 1:400]
#dim(my.data.mini)
# ### CHOOSE DATASET For testing
# data <- my.data.mini #mini
data <- my.data.hfs

# Calculate genetic distance using Weir & Cockerham Fst
genet.dist(dat = data, method = "WC84")

# Provides Hobs (Ho), mean gene diversities w/in pops (Hs), Fis, Fst
basicstat <- basic.stats(data, diploid=T, digits=2)
names(basicstat) 

# lots of great data and info here
summary(basicstat$perloc[, "Fst"]) # doesn't really work

plot(basicstat$perloc[,"Fst"]) # plot Fst by index

wc(data)
# Fst following Weir and Cockerham's estimate


# which markers are greater than a specific Fst value?
val <- 0.1
rownames(basicstat$perloc)[which(basicstat$perloc[, "Fst"] > val)]
high.fst.markers <- rownames(basicstat$perloc)[which(basicstat$perloc[, "Fst"] > val)]
length(high.fst.markers)

# extract these from the genlight object and re-run nj
head(locNames(my.data.gid))
high.fst.markers.true <- gsub(x = high.fst.markers, pattern = "^X", replacement = "", perl = T)

my.data.high.fst <- my.data.gid[, loc=high.fst.markers.true]
my.data.high.fst


# Try dapc now
dapc1 <- dapc(my.data.high.fst, n.pca = 10, n.da = 1) 
# n.pca = number axes to be retained; n.da = number of axes retained in Discriminant Analysis step
dapc1 # variance = 0.3

# Density plot of samples along discriminant function 1
scatter(dapc1, scree.da = F, bg = "white", posi.pca = "topright", legend = T
        #, txt.leg=c("Pipestem","Pendrell","DeepBay","Rosewall","BayneSpat","PendrellFarm")
)

# Composition plot
compoplot(dapc1
          #, lab="" # comment out if you want sample labels
          , txt.leg=c("Pipestem","Pendrell","DeepBay","Rosewall","BayneSpat","PendrellFarm")
)

# Loading plot
loadingplot(dapc1$var.contr, thres=1e-3)




#### TODO ####
# 1. Use seppop to separate populations then compare (e.g. pendrell vs pendrell farm; pipestem vs pendrell)
# 2. How to get statistics on Fst?
# 3. How to calculate pi across genome?
# 4. How to best incorporate aligned markers against the Pacific Oyster chromosome genome assembly here? 
# 5. How to best calculate population-specific metrics, should I just use Stacks output?
















### DEAD CODE ###

# seploc returns a list of objects, each corresp. to a marker

# sep.gid <- seploc(my.data.gid)
# class(sep.gid)
# head(names(sep.gid))
# 
# # can repool this way:
# new.obj <- NULL
# for(i in 1:length(high.fst.markers)){
#   high.fst.markers[i]
# }
# 
# repool()
# 
# 

# ### NOT WORKING ###
# # Hierarchical Fst tests (AMOVA for SNP datasets)
# loci <- data[,-1] # remove the population column
# populations <- as.character(data[,1]) # select only the population column
# 
# ### The following use hierarchies, so we need more than one distinguisher
# # not working: varcomp.glob(levels = data.frame(populations), loci = loci, diploid = T)
# 
# # not working: varcomp.glob(loci=data) # return multilocus estimators of variance components and Fstats
# 
# # not working: test.g(loci, level = populations)
#