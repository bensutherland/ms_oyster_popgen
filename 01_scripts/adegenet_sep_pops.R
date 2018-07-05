# Separate Analysis of individual populations

#### Front Matter ####
# Clean space
# rm(list=ls())

# Set working directory
# Xavier
setwd("~/Documents/01_moore_oyster_project/stacks_workflow_no_Qc")

# Load part 1 results
load("11-other_stats/adegenet_output.RData")

my.data.gid

# Separate gid
sep.obj <- seppop(x = my.data.gid)
names(sep.obj)

wild.bc.gid <- repool(sep.obj$Hisnit, sep.obj$Pendrell, sep.obj$PendrellFarm, sep.obj$Pipestem, sep.obj$Serpentine)
bc.gid <- repool(sep.obj$Pendrell, sep.obj$Hisnit)
fr.gid <- repool(sep.obj$FranceW, sep.obj$FranceC)
ch.gid <- repool(sep.obj$ChinaQDW, sep.obj$ChinaBe) 
# ChinaQDW = spat taken from culture region (QD) by artificial attachment, then grown in 'lantern net culture'
# ChinaBe = from QD region, sampled directly from rock on the beach

#### Choose dataset ####
data.gid <- wild.bc.gid
# data.gid <- bc.gid
# data.gid <- fr.gid
# data.gid <- ch.gid

#### 1. Investigate Wild Populations Only ####
# Identify the populations in the data
levels(pop(x = data.gid))
nPop(data.gid)
nInd(data.gid)
indNames(data.gid)

# Show sample size per population
par(mfrow=c(1,1), mar=c(8,5,3,3))
barplot(table(pop(data.gid)), col=funky(17)
        , las=2
        , xlab=""
        , ylab="Sample size"
        , ylim = c(0,40))
abline(h = c(10,20,30), lty=2)

# Change from genind file to hfstat
data.hf <- genind2hierfstat(data.gid)
rownames(data.hf) <- indNames(data.gid)

# PCA on a matrix of individual genotype frequencies (hierfstat)
y <- indpca(data.hf, ind.labels = rownames(data.hf))
plot(y, cex = 0.7)

# DAPC
dapc <- dapc(data.gid, n.pca = 10, n.da = 1)
scatter(dapc, scree.da = F, bg = "white", legend = T
        , txt.leg=rownames(dapc$means)
        , posi.leg = "topleft"
)

# Composition plot (barplot showing the probabilities of assignments of individuals to the different clusters)
par(mar=c(10,3,3,3))
compoplot(dapc
          , txt.leg = rownames(dapc$means)
          , posi = "topright"
          , cex.names = 0.6 
)

# Loading plot # Plot marker variance contribution to DAPC
par(mfrow=c(1,1), mar=c(3,4,3,3))

# Plot the loading values of the different markers into the DAPC
loadingplot(dapc$var.contr, thres=1e-3, las = 1)

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
pairwise.WCfst(data.hf)

# Bootstrapping
# requires latest hierfstat (v0.04-29) otherwise get error
library(devtools)
install_github("jgx65/hierfstat")
library("hierfstat")
boot.fst.wild.bc <- boot.ppfst(dat = data.hf, nboot = 100, quant = c(0.025,0.975))
boot.fst.wild.bc


#### 8. BC Seq Facility (Hisnit, Pendrell) ####

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
