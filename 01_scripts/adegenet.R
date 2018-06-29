# Import data from PLINK, conduct basic pop gen analyses, create a neighbour-joining tree, PCA, DAPC

#### Front Matter ####
# Clean space
# rm(list=ls())

# Set working directory
# Xavier
setwd("~/Documents/01_moore_oyster_project/stacks_workflow_no_Qc")

# Wayne
# setwd("~/Documents/miller/00_Moore_oyster_project/04_analysis/stacks_workflow")

# Resource packages required:
#install.packages("purrr")
library("purrr")
#install.packages("dplyr")
library("dplyr")

# require(devtools)
# install_version("hierfstat", version = "0.04-22", repos = "http://cran.us.r-project.org") # compoplot functional
# install_github("jgx65/hierfstat") # compoplot not functional
library(hierfstat)

library("ape")
library("pegas")
library("seqinr")
library("ggplot2")

# install_version("adegenet", version = "2.0.1", repos = "http://cran.us.r-project.org") # compoplot functional
# install.packages("adegenet") # compoplot not functional
library("adegenet")

sessionInfo()

#### 1. Import Genomics (genlight) ####
my.data <- read.PLINK(file = "06-stacks_rx/batch_1.raw")
# my.data <- read.PLINK(file = "07-filtered_vcfs/batch_1_0.5Mreads_r70_p7_maf0.05.raw")
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
# Apply here if needed (e.g. remove samples, populations etc.)

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
par(mfrow=c(1,1), mar=c(3,3,3,3))
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
# or, if missing data:
# plot(njs(D), type ="fan", cex = 0.7)

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
text(x = -2, y = 11, paste(labels=nLoc(my.data), "loci", sep = " "))
text(x = -2, y = 9, labels=paste("PC1=Red", "\n", "PC2=Green", "\n", "PC3=Blue"))

# Bring colorplot colors into the samples of the Neighbour-joining tree
plot(nj(D), type="fan", show.tip=T, cex =  0.75)
# or
# plot(njs(D), type="fan", show.tip=T, cex =  0.8)
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
dapc1 # variance = 0.112

# Density plot of samples along discriminant function 1
scatter(dapc1, scree.da = F, bg = "white", legend = T
        , txt.leg=rownames(dapc1$means)
        , posi.leg = "topleft"
        )

# assignplot hasn't been explored/implemented here yet
#assignplot(x = dapc1)

# Composition plot (barplot showing the probabilities of assignments of individuals to the different clusters)
par(mar=c(10,3,3,3))
#pdf(file = "test", width = 10, height = 10)
compoplot(dapc1
          #, lab="" # comment out if you want sample labels
          , txt.leg = rownames(dapc1$means)
          , posi = "topright"
#          , cex = 0.7
          , cex.names = 0.6 
          )
#dev.off()

# Loading plot # Plot marker variance contribution to DAPC
par(mfrow=c(1,1), mar=c(3,4,3,3))

# Plot the loading values of the different markers into the DAPC
loadingplot(dapc1$var.contr, thres=1e-3)

# Save out the loading values for the dapc linked to the marker names
par(mar=c(5,5,4,4))
dapc1$pca.loadings[1:5,]
max(dapc1$pca.loadings[,1])
max(dapc1$var.contr) # this is what's being plotted
plot(dapc1$pca.loadings[,1], dapc1$var.contr ) # absolute value of pca loadings are 
plot(abs(dapc1$pca.loadings[,1]), dapc1$var.contr ) # absolute value of pca loadings are highly correlated to the 'var.contrib'
# Where are the locus names? They aren't included in this, which seems to just take the 'index' of the locus. 
# Therefore, we have to currently assume that the order is the same as in my.data
locus.and.dapc1.var.contrib <- as.data.frame(cbind(locNames(my.data), dapc1$var.contr))
colnames(locus.and.dapc1.var.contrib) <- c("mname","dapc1.var.contr")
write.csv(x = locus.and.dapc1.var.contrib, file = "11-other_stats/locus_and_dapc1_var_contrib.csv"
          , quote = F, row.names = F)


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
par(mfrow=c(1,1), mar=c(8,5,3,3))
barplot(table(pop(my.data.gid)), col=funky(17)
        #, las=3, las = 1
        , las=2
        , xlab=""
        , ylab="Sample size"
        , ylim = c(0,40))
abline(h = c(10,20,30), lty=2)

# How many alleles? (should all be 2)
table(nAll(my.data.gid))

# Change from genind to hierfstat
all.data.hf <- genind2hierfstat(my.data.gid)
rownames(all.data.hf) <- indNames(my.data.gid)

# Pairwise Fst
pairwise.WCfst(all.data.hf)

# Bootstrapping
# requires latest hierfstat (v0.04-29) otherwise get error
# library(devtools)
# install_github("jgx65/hierfstat")
# library("hierfstat")
boot.fst.all <- boot.ppfst(dat = all.data.hf, nboot = 100, quant = c(0.025,0.975))
boot.fst.all

write.csv(boot.fst.all$ll, file = "boot_fst_all_ll.csv")
write.csv(boot.fst.all$ul, file = "boot_fst_all_ul.csv")

#### 7. Investigate Wild Populations Only ####
sep.obj <- seppop(x = my.data.gid)
names(sep.obj)

# Wild populations (and moved to farm)
wild.bc.gid <- repool(sep.obj$Hisnit, sep.obj$Pendrell, sep.obj$PendrellFarm, sep.obj$Pipestem, sep.obj$Serpentine)
nPop(wild.bc.gid)
nInd(wild.bc.gid)
indNames(wild.bc.gid)

# Show sample size per population
par(mfrow=c(1,1), mar=c(8,5,3,3))
barplot(table(pop(wild.bc.gid)), col=funky(17)
        #, las=3, las = 1
        , las=2
        , xlab=""
        , ylab="Sample size"
        , ylim = c(0,40))
abline(h = c(10,20,30), lty=2)

# Change from genind file to hfstat
wild.bc.hf <- genind2hierfstat(wild.bc.gid)
rownames(wild.bc.hf) <- indNames(wild.bc.gid)

# Perform PCA on a matrix of individual genotype frequencies (hierfstat)
y <- indpca(wild.bc.hf, ind.labels = rownames(wild.bc.hf))
plot(y, cex = 0.7)

# Perform dapc on the wild only individuals
dapc2 <- dapc(wild.bc.gid, n.pca = 10, n.da = 1)
scatter(dapc2, scree.da = F, bg = "white", legend = T
        , txt.leg=rownames(dapc2$means)
        , posi.leg = "topleft"
        )

# Composition plot (barplot showing the probabilities of assignments of individuals to the different clusters)
par(mar=c(10,3,3,3))
compoplot(dapc2
          #, lab="" # comment out if you want sample labels
          , txt.leg = rownames(dapc2$means)
          , posi = "topright"
          #          , cex = 0.7
          , cex.names = 0.6 
)

# Loading plot # Plot marker variance contribution to DAPC
par(mfrow=c(1,1), mar=c(3,4,3,3))

# Plot the loading values of the different markers into the DAPC
loadingplot(dapc2$var.contr, thres=1e-3)

# Identify the top loading marker names
toploading.markers <- dimnames(dapc2$var.contr)[[1]][which(dapc2$var.contr > 0.0005)]
toploading.markers <- gsub(pattern = "\\..*", replacement = "", toploading.markers) # match the genlight name format (remove .1 or .2)
toploading.markers

# # subset the gid file w/ toploading markers
# locNames(wild.bc.gid)
# 
# wild.bc.gid.toploading <- wild.bc.gid[, loc=toploading.markers]
# wild.bc.gid.toploading


#### TOPLOADING STUFF ####
# Subset the original my.data genlight file with the toploading markers and convert back to a genlight file
my.data.mat <- as.matrix(my.data)
my.data.mat[1:5,1:5]
# retain only top loading markers
my.data.mat.toploading <- my.data.mat[,toploading.markers]
# retain only wild samples (removing DeepBay only currently)
samples.to.keep <- rownames(my.data.mat)[grep(pattern = "DeepBay", x = rownames(my.data.mat), invert = T)]
my.data.mat.toploading.wild <- my.data.mat.toploading[samples.to.keep,]

# Convert back to genlight format
my.data.toploading.gl  <- as.genlight(my.data.mat.toploading.wild)

# Determine distance between individuals
# First use seploc to create a list of ten blocks of loci, each to be computed in parallel
x <- seploc(my.data.toploading.gl, n.block=10, parallel=F)

## Use dist in a lapply loop to compute pairwise distances between individuals for each block
lD <- lapply(x, function(e) dist(as.matrix(e)))

# Obtain general distance matrix by summing these
D <- Reduce("+", lD)

# Plot a neighbor-joining tree using ape's nj function on the distance matrix
plot(nj(D), type="fan", cex=0.7)
# this fails if there are too few markers
# plot(njs(D), type="fan", cex=0.7)

# Pairwise Fst
pairwise.WCfst(wild.bc.hf)

# Bootstrapping
# requires latest hierfstat (v0.04-29) otherwise get error
library(devtools)
install_github("jgx65/hierfstat")
library("hierfstat")
boot.fst.wild.bc <- boot.ppfst(dat = wild.bc.hf, nboot = 100, quant = c(0.025,0.975))
boot.fst.wild.bc


#### 8. BC Seq Facility (Hisnit, Pendrell) ####
bc.gid <- repool(sep.obj$Pendrell, sep.obj$Hisnit)
# could also use # sep.obj$Serpentine
nPop(bc.gid)
nInd(bc.gid)
indNames(bc.gid)
pop(bc.gid)

# If you want to edit the population variable, do it here
indiv.used <- as.numeric(gsub(pattern = "*.*\\_", replacement = "", indNames(bc.gid)))
indiv.used <- as.data.frame(indiv.used)
indNames(bc.gid)
refined.pop <- c(rep("Pen_2", times = 15), rep("Pen_3", times = 13), rep("His_3", times = 5), rep("His_6", times = 11)
                 , rep("His_3", times = 5))
# TODO # import sample database for auto matching of this using a merge function w/ sample ID and size class

length(indNames(bc.gid))
length(refined.pop)

pop(bc.gid) <- refined.pop
pop(bc.gid)

# DAPC
dapc3 <- dapc(bc.gid, n.pca = 10, n.da = 1)
scatter(dapc3, scree.da = F, bg = "white", legend = T
        , txt.leg=rownames(dapc3$means)
        , posi.leg = "topleft"
)
# var (proportion of conserved variance): 0.244


# Composition plot (barplot showing the probabilities of assignments of individuals to the different clusters)
par(mar=c(10,3,3,3))
compoplot(dapc3
          #, lab="" # comment out if you want sample labels
          , txt.leg = rownames(dapc3$means)
          , posi = "topright"
          #          , cex = 0.7
          , cex.names = 0.6 
)

# Loading plot # Plot marker variance contribution to DAPC
par(mfrow=c(1,1), mar=c(3,4,3,3))

# Plot the loading values of the different markers into the DAPC
loadingplot(dapc3$var.contr, thres=1e-3)

# Show sample size per population
par(mfrow=c(1,1), mar=c(8,5,3,3))
barplot(table(pop(bc.gid)), col=funky(17)
        #, las=3, las = 1
        , las=2
        , xlab=""
        , ylab="Sample size"
        , ylim = c(0,40))
abline(h = c(10,20,30), lty=2)


# Change to hierfstat
bc.hf <- genind2hierfstat(bc.gid)
rownames(bc.hf) <- indNames(bc.gid)

# PCA
z <- indpca(bc.hf, ind.labels = rownames(bc.hf))
# z <- indpca(bc.hf, ind.labels = rownames(refined.pop)) # use this if want to see the refined pop
# plot PCA
plot(z, cex = 0.7)


# Basic stats (gives per loc and overall, but only one value, not one for each comparison)
bc.stats <- basic.stats(bc.hf, diploid = T, digits = 4)

par(mar=c(5,5,2.5,3))
plot(bc.stats$perloc$Fst, ylab = "Fst", las = 1) # Plot Fst
text(x = 2000, y = 0.2, labels = paste(nLoc(my.data),"loci" ))
text(x = 2000, y = 0.19, labels = paste(nInd(my.data),"ind" ))
text(x = 2000, y = 0.18, labels = "Size class, Hisnit-Serp")
text(x = 2000, y = 0.17, labels = "Max. Fst")


# Save out these Fst values per locus to match to location in the genome
perloc.stats.df <- bc.stats$perloc # save df w/ Ho, Hs, Ht, Dst, Htp, Dstp, Fst, Fstp, Fis, Dest (from hierfstat)
write.csv(x = perloc.stats.df, file = "../ms_oyster_popgen/perloc.stats.csv")
# match these to their genomic location and plot in Rscript (#todo)

# Weir-Cockerham
pairwise.wcfst.bc.hf <- pairwise.WCfst(dat = bc.hf)

# Bootstrapping
boot.fst.obj <- boot.ppfst(dat = bc.hf, nboot = 1000, quant = c(0.025,0.975))
boot.fst.obj


#### 9. Hobs vs Hexp ####
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


### HERE BE DRAGONS ###

# The following can be used to test for significant difference between Hobs and Hexp
# # test : H0: Hexp = Hobs (using Bartlett test of homogeneity of variances)
# exp.v.obs.test <- bartlett.test(list(div$Hexp, div$Hobs))
# exp.v.obs.test
# options(scipen=3)
# text(x = 0.6, y = 0.9, labels = paste("pval =", exp.v.obs.test$p.value, sep = ""))


# #### 10. Hardy-Weinberg Equilibrium ####
# # Choose the dataset to test
# #dataset.for.hwt <- my.data.gid
# dataset.for.hwt <- wild.bc.nofarm.gid
# 
# # run the test
# my.data.hwt <- hw.test(dataset.for.hwt, B=0) 
# head(my.data.hwt)
# table(my.data.hwt[,3] < 0.01) 
# # for all samples (95) incl. farms  = False 4906; True 6364
# # for wild.bc.nofarm.gid            = False 6270; True 4998
# # Note, this looks as though it is being tested regardless of population, and remember
# # When there is population structure, there is often loci out of hwt
# # However, even when looking at samples without pop structure, it appears that many loci are out of HWE



##### 11. Other Fstats ####
# May also try diveRsity or StAMPP

# Another method, conduct pairwise fst on the data (using Nei's estimator)
my.data.fst <- pairwise.fst(x = my.data.gid, pop = pop(my.data.gid))
my.data.fst
is.euclid(my.data.fst) # False, there are zero distances (probably the negative values!)

# Provide three F statistics (Fst (pop/total), Fit (Ind/total), Fis (ind/pop))
# fstat(my.data.gid)

# which markers are greater than a specific Fst value?
val <- 0.1
rownames(basicstat$perloc)[which(basicstat$perloc[, "Fst"] > val)]
high.fst.markers <- rownames(basicstat$perloc)[which(basicstat$perloc[, "Fst"] > val)]
length(high.fst.markers)

# extract these from the genlight object
head(locNames(my.data.gid))
high.fst.markers.true <- gsub(x = high.fst.markers, pattern = "^X", replacement = "", perl = T)

my.data.high.fst <- my.data.gid[, loc=high.fst.markers.true]
my.data.high.fst


