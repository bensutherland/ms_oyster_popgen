# This script is used to import data from plink, conduct basic analyses, create a neighbour-joining tree, PCA, DAPC
# then 

#### Front Matter ####
# Clean space
# rm(list=ls())

# Set working directory
# Xavier
setwd("~/Documents/01_moore_oyster_project/stacks_workflow_2018-05-03/")

# Wayne
# setwd("~/Documents/miller/00_Moore_oyster_project/04_analysis/stacks_workflow")

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

#### 1. Import Genomics (genlight) ####
my.data <- read.PLINK(file = "06-stacks_rx/batch_1.raw")
my.data

## Explore data attributes
nInd(my.data) # how many individuals?
nLoc(my.data) # how many loci?
indNames(my.data) # see names of individuals
nPop(my.data) # how many pops?
pop(my.data) # what are the pop names?
head(locNames(my.data, withAlleles = T), n = 20) # see allele names

# Rename pop attributes
pop(my.data) <- gsub(x = indNames(my.data), pattern = "\\_.*", replacement = "") # create pop ID
pop(my.data) # what are the pop names?

#### 1.b. Additional filtering ####
# Todo: how to remove populations from the analysis?

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
par(mfrow=c(1,1), mar=c(3,3,2,3))
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
plot(nj(D), type="fan", cex=0.7)


#### 4. Principle Components Analysis ####
pca1 <- glPca(my.data, nf = 3)

# Plot PCA (axes 1,2,3)
par(mfrow=c(2,1))
scatter(x = pca1, posi = "topleft", xax = 1, yax = 2)
title("PC1 (x) vs PC2 (y)", adj = 1)

scatter(x = pca1, posi = "topleft", xax = 3, yax = 2) # show PC3 and PC4
title("PC3 (x) vs PC2 (y)", adj = 1)

# Plot allele loadings
num.retained.pcs <- length(dimnames(pca1$loadings)[[2]]) # how many axes were retained? 

# Plot loading values of the markers in the PCs
par(mfrow=c(num.retained.pcs,1))
# Plot the loading values of the different markers into the PCA
for(i in 1:num.retained.pcs){
  loadingplot(x = pca1, axis = i
              , main = paste("PC",i,"_loadings_(alleles)", sep = "")
              #, threshold = quantile(x = pca1$loadings[,i], 0.8) # not working
              , threshold = 0.001
              )
}

# Plot colorplot using the PC scores of the samples
par(mfrow=c(1,1))
myCol <- colorplot(pca1$scores, pca1$scores, transp=T, cex=4)
abline(h=0,v=0, col = "grey")
text(x = -2, y = 15, paste(labels=nLoc(my.data), "loci", sep = " "))
text(x = -2, y = 12, labels=paste("PC1=Red", "\n", "PC2=Green", "\n", "PC3=Blue"))

# Bring colorplot colors into the samples of the Neighbour-joining tree
plot(nj(D), type="fan", show.tip=T, cex =  0.8)
tiplabels(pch=20, col=myCol, cex=4)


# ## try w/ smaller subset of samples
# pop(my.data)
# sep.obj <- seppop(x = my.data)
# names(sep.obj)
# pendrell.pipestem.gid <- repool(sep.obj$Pendrell, sep.obj$Pipestem) # make a pendrell pipestem only dataset
# nPop(pendrell.pipestem.gid)
# nInd(pendrell.pipestem.gid)



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
par(mar=c(10,3,3,3))
compoplot(dapc1
          #, lab="" # comment out if you want sample labels
          , txt.leg=c("Pipestem","Pendrell","DeepBay","Rosewall","BayneSpat","PendrellFarm")
          , posi = "topright"
          #, cex.lab = 0.7 #doesn't seem to work unfortunately
          )

# Loading plot # Plot marker variance contribution to DAPC
par(mfrow=c(1,1), mar=c(3,4,3,3))

# Plot the loading values of the different markers into the PCA
loadingplot(dapc1$var.contr, thres=1e-3)


####6. Convert genlight to genind object (indiv. genotypes) ####
# First, convert to matrix
my.data.mat <- as.matrix(my.data)
my.data.mat[1:5,1:5]

# Second, translate the number of minor allele to genind format
my.data.mat[my.data.mat == 0] <- "1/1" #homozygote reference
my.data.mat[my.data.mat == 1] <- "1/2" #heterozygote
my.data.mat[my.data.mat == 2] <- "2/2" #homozygote alternate
#my.data.mat[my.data.mat == "NA"] <- "NA" # NA
my.data.mat[1:5,1:5]

# Convert matrix to genind
my.data.gid <- df2genind(my.data.mat, sep = "/", ploidy = 2) # convert df to genind

# Generate populations for genind
pop(my.data.gid) # currently null
indNames(my.data.gid) # The individual names actually contain pop IDs
pop(my.data.gid) <- gsub(x = indNames(my.data.gid), pattern = "\\_.*", replacement = "") # create pop ID
pop(my.data.gid)

# Show sample size per population
table(pop(my.data.gid))
par(mfrow=c(1,1), mar=c(6,5,3,3))
barplot(table(pop(my.data.gid)), col=funky(17), las=3, las = 1
        ,  xlab="Population", ylab="Sample size"
        , ylim = c(0,40))

# How many alleles? (should all be 2)
table(nAll(my.data.gid))


##### Separate Populations and do some PCA #####
sep.obj <- seppop(x = my.data.gid)
names(sep.obj)

# WILD BC AND MOVED TO FARM
wild.bc.gid <- repool(sep.obj$BayneSpat, sep.obj$Pendrell, sep.obj$PendrellFarm, sep.obj$Pipestem)
nPop(wild.bc.gid)
nInd(wild.bc.gid)
indNames(wild.bc.gid)
#change to hf
wild.bc.hf <- genind2hierfstat(wild.bc.gid)
rownames(wild.bc.hf) <- indNames(wild.bc.gid)

y <- indpca(wild.bc.hf, ind.labels = rownames(wild.bc.hf))
plot(y, cex = 0.7)

# WILD BC NO FARM
wild.bc.nofarm.gid <- repool(sep.obj$BayneSpat, sep.obj$Pendrell, sep.obj$Pipestem)
nPop(wild.bc.nofarm.gid)
nInd(wild.bc.nofarm.gid)
indNames(wild.bc.nofarm.gid)
#change to hf
wild.bc.nofarm.hf <- genind2hierfstat(wild.bc.nofarm.gid)
rownames(wild.bc.nofarm.hf) <- indNames(wild.bc.nofarm.gid)

z <- indpca(wild.bc.nofarm.hf, ind.labels = rownames(wild.bc.nofarm.hf))
plot(z, cex = 0.7)






## PENDRELL AND PIPESTEM ONLY
pendrell.pipestem.gid <- repool(sep.obj$Pendrell, sep.obj$Pipestem) # make a pendrell pipestem only dataset
# note: expect warning message that non-type markers are deleted
nPop(pendrell.pipestem.gid)
nInd(pendrell.pipestem.gid)

basicstat.pendrell.pipestem <- basic.stats(pendrell.pipestem.gid, diploid = T, digits = 4)
names(basicstat.pendrell.pipestem) #n indiv, pop freq, Ho, Hs, Fis, perloc, overall
str(basicstat.pendrell.pipestem)

plot(basicstat.pendrell.pipestem$perloc$Fst)

# #bootstrapping NOT WORKING
# bootstrap.pipestem.pendrell <- boot.ppfst(pendrell.pipestem.gid)

indNames(pendrell.pipestem.gid)
pendrell.pipestem.hf <- genind2hierfstat(pendrell.pipestem.gid)
dim(pendrell.pipestem.hf)
rownames(pendrell.pipestem.hf) <- indNames(pendrell.pipestem.gid)

x <- indpca(pendrell.pipestem.hf, ind.labels = rownames(pendrell.pipestem.hf))
plot(x, cex = 0.7)




#










#### 7. Hobs vs Hexp ####

# Descriptive Statistics on genind file
div <- summary(my.data.gid) # genetic diversity w/ adegenet
str(div) # What is contained within this summary object
# contains: sample size, sample size by pop, number loci and names, population IDs, % NA, Hobs, Hexp

# Plotting Hobs by Hexp
par(mfrow=c(1,3), mar=c(4,3,3,3))
plot(div$Hobs
     , xlab="Locus index", ylab="Obs. Het."
     , main = "Observed heterozygosity per locus"
     , ylim = c(0,1))
# Note that stacks calculates obs het diff than other programs.. (missing data is not considered properly)
plot(div$Hexp
     , xlab="Locus index", ylab="Exp. Het."
     , main = "Expected heterozygosity per locus"
     , ylim = c(0,1))

# Plot as a function of each other:
plot(div$Hobs, div$Hexp
     , xlab = "Obs. Het.", ylab="Exp. Het."
     , main = "Expected Heterozygosity vs. Observed"
     , ylim = c(0,1))
abline(0, 1, lty=2)


# # Plot inverse to above
# plot(div$Hexp, div$Hobs
#      , pch=20, cex=3, xlim=c(0, 0.6), ylim=c(0, 1)
#      , main = paste("nLoc: ", nLoc(my.data.gid),"nInd: "
#                     , nInd(my.data.gid), "nPop: ", nPop(my.data.gid)))
# abline(0,1,lty=2)

# The following can be used to test for significant difference between Hobs and Hexp
# # test : H0: Hexp = Hobs (using Bartlett test of homogeneity of variances)
# exp.v.obs.test <- bartlett.test(list(div$Hexp, div$Hobs))
# exp.v.obs.test
# options(scipen=3)
# text(x = 0.6, y = 0.9, labels = paste("pval =", exp.v.obs.test$p.value, sep = ""))


#### Test for significance on genind file ####

#### Hardy-Weinberg Equilibrium ####
# Choose the dataset to test
#dataset.for.hwt <- my.data.gid
dataset.for.hwt <- wild.bc.nofarm.gid

# run the test
my.data.hwt <- hw.test(dataset.for.hwt, B=0) 
head(my.data.hwt)
table(my.data.hwt[,3] < 0.01) 
# for all samples (95) incl. farms  = False 4906; True 6364
# for wild.bc.nofarm.gid            = False 6270; True 4998
# Note, this looks as though it is being tested regardless of population, and remember
# When there is population structure, there is often loci out of hwt
# However, even when looking at samples without pop structure, it appears that many loci are out of HWE



##### 7. Fstats ####
# May also try diveRsity or StAMPP

# Conduct pairwise fst on the data (using Nei's estimator)
my.data.fst <- pairwise.fst(x = my.data.gid, pop = pop(my.data.gid))
my.data.fst
is.euclid(my.data.fst) # Flase, there are zero distances (probably the negative values!)


# Provde three F statistics (Fst (pop/total), Fit (Ind/total), Fis (ind/pop))
fstat(my.data.gid)








# Create a hierfstat object (nloci+1 = columns; ninds = rows)
my.data.hfs <- genind2hierfstat(my.data.gid)
class(my.data.hfs) # is a dataframe
dim(my.data.hfs) # with 95 indiv and 11,271 loci

# Name the rows with the sample IDs
rownames(my.data.hfs) <- indNames(my.data.gid)

# OPTIONAL TO LIMIT DATA
#my.data.mini <- my.data.hfs[, 1:400]
#dim(my.data.mini)
# ### CHOOSE DATASET For testing
# data <- my.data.mini #mini

### Working w/ hierfstat obj
data <- my.data.hfs

# define populations for this data
pops <- gsub(x = rownames(data), pattern = "\\_.*", replacement = "", perl = T)







## OLDER ##
# Calculate genetic distance using Weir & Cockerham Fst
genet.dist(dat = data, method = "WC84")

# Test signif of fx of test.lev on genetic diffn
test.between(data = data, test.lev = pops, nperm = 2, rand.unit = )

# Test the significance of the effect of level on genetic differentiation
test.g(data = data, level = pops, nperm = 2)

# 
my.obj <- boot.ppfst(dat = data, nboot = 2)




# Provides Hobs (Ho), mean gene diversities w/in pops (Hs), Fis, Fst
basicstat <- basic.stats(data, diploid=T, digits=2)
names(basicstat) 

# lots of great data and info here
summary(basicstat$perloc[, "Fst"]) # doesn't really work

plot(basicstat$perloc[,"Fst"]) # plot Fst by index

# wc(data)
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