# Import data from PLINK, conduct basic pop gen analyses, create a neighbour-joining tree, PCA, DAPC

#### Front Matter ####
# Clean space
# rm(list=ls())

## Install and load packages
# install.packages("purrr")
# install.packages("dplyr")
# install.packages("ape")
# install.packages("pegas")
# install.packages("seqinr")
# install.packages("ggplot2")
# install.packages("devtools")

# Specifics for hierfstat and adegenet
# require(devtools)
# install_version("hierfstat", version = "0.04-22", repos = "http://cran.us.r-project.org") # compoplot functional
# install_github("jgx65/hierfstat") # compoplot not functional

# install.packages("hierfstat")
library("hierfstat")

# install_version("adegenet", version = "2.0.1", repos = "http://cran.us.r-project.org") # compoplot functional
# install.packages("adegenet") # compoplot not functional
library("adegenet")


library("purrr")
library("dplyr")
library("ape")
library("pegas")
library("seqinr")
library("ggplot2")
library("devtools")


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


#### 1. Import data ####
my.data <- read.PLINK(file = "11-adegenet_analysis/populations_single_snp_HWE.raw")
my.data

## Explore data attributes
nInd(my.data) # How many individuals?
nLoc(my.data) # How many loci?
indNames(my.data)[1:20] # What individuals?
nPop(my.data) # How many pops?
unique(pop(my.data)) # What pops?
locNames(my.data, withAlleles = T)[1:20] # What allele names?

#### 1.1 Edit data attributes
## Rename populations
# example to edit populations
# pop(my.data) <- gsub(x = indNames(my.data), pattern = "\\_.*", replacement = "") # create pop ID


#### 2. Basic Analysis (adegenet) ####
# Plot instances of minor allele across individuals and loci
png(file = paste0(output.dir, "glPlot_all.png"), width = 924, height = 600)
glPlot(x = my.data, posi="topleft") 
dev.off()

# Create density plot of minor allele frequencies
myFreq <- glMean(my.data)
pdf(file = paste0(output.dir, "maf_hist_all.pdf"), width = 6, height = 4)
hist(myFreq, proba=T, col="gold", xlab = "Allele frequencies"
     , main = ""
     , ylim = c(0,50)
     , ylab = "Density of second allele frequencies"
     )
text(x = 0.4, y = 7, labels = paste(nLoc(my.data), " loci", sep = "" ))
temp <- density(myFreq)
lines(temp$x, temp$y, lwd=3)
dev.off()


#### 3. Neighbour-joining tree, all loci ####
par(mfrow=c(1,1), mar=c(3,3,3,3))

## Calculate distance between individuals
# Create blocks of loci for parallelization
x <- seploc(my.data, n.block=10, parallel=F)

## Use dist in a lapply loop to compute pairwise distances between individuals for each block
lD <- lapply(x, function(e) dist(as.matrix(e)))
class(lD) # list
names(lD) # 10 blocks
class(lD[[1]]) # lD is a list of length 10 containing distance matrices bw pairs of indiv

# Obtain general distance matrix by summing these
D <- Reduce("+", lD)

## Plot distance
# Neighbor-joining tree (NJT) with ape
pdf(file = paste0(output.dir, "nj_tree_all.pdf"), width = 11, height = 9)
plot(nj(D), type="fan", cex=0.7)
dev.off()

# In the case of missing data
# plot(njs(D), type ="fan", cex = 0.7)


#### 4. Principle Components Analysis, all loci ####
# Perform PCA
pca1 <- glPca(my.data, nf = 3)

# Plot PCA (axes 1,2,3)
pdf(file = paste0(output.dir, "pca_all_samples.pdf"), width = 11.5, height = 7.5)
par(mfrow=c(2,1))
scatter(x = pca1, posi = "topleft", xax = 1, yax = 2)
title("PC1 (x) vs PC2 (y)", adj = 1)

scatter(x = pca1, posi = "topleft", xax = 3, yax = 2) # show PC3 and PC4
title("PC3 (x) vs PC2 (y)", adj = 1)
dev.off()

# Plot allele loadings
num.retained.pcs <- length(dimnames(pca1$loadings)[[2]]) # how many axes were retained? 

# Plot loading values of the markers in the PCs
pdf(file = paste0(output.dir, "pc_loadings.pdf"), width = 8, height = 8)
par(mfrow=c(num.retained.pcs,1))
# Plot the loading values of the different markers into the PCA
for(i in 1:num.retained.pcs){
  loadingplot(x = pca1, axis = i
              , main = paste("PC",i,"_loadings_(alleles)", sep = "")
              #, threshold = quantile(x = pca1$loadings[,i], 0.8) # not working
              , threshold = 0.001
              )
}
dev.off()

# Plot PCA with samples coloured by PC
pdf(file = paste0(output.dir , "pca_colorplot.pdf"), width = 8, height = 8)
par(mfrow=c(1,1), mar = c(4,4,4,4))
myCol <- colorplot(pca1$scores, pca1$scores, transp=T, cex=4)
abline(h=0,v=0, col = "grey")
text(x = 2, y = 5, paste(labels=nLoc(my.data), "loci", sep = " "))
text(x = 2, y = 4, labels=paste("PC1=Red", "\n", "PC2=Green", "\n", "PC3=Blue"))
dev.off()

# Bring colorplot colors into the samples of the Neighbour-joining tree
pdf(file = paste0(output.dir, "njt_colorplot.pdf"), width = 11, height = 9)
plot(nj(D), type="fan", show.tip=T, cex =  0.75)
tiplabels(pch=20, col=myCol, cex=4)
dev.off()


#### to remove ####
# Define colours
# my.pops <- levels(dapc1$grp)
# my.cols <- c("darkseagreen4", "darkseagreen3" #china1
#              , "blue3" #dpb
#              , "gold1", "gold3" #FRA and FRAF (darker = derived)
#              , "darkslategray2" #GUR
#              , "mediumpurple3" #HIS
#              , "turquoise4" #JPN
#              , "purple1", "purple4" #PEN PENF
#              , "darkviolet" #PIP 
#              ,"darkseagreen1" #QDC
#              , "red" #ROS
#              , "darkseagreen" #china2
#              , "orchid2"
#              )
# 
# my.cols.df <- as.data.frame(cbind(my.pops, my.cols))
# write.csv(my.cols.df, file = "01-info_files/my_cols.csv", quote = F, row.names = F)

#### 5. Discriminant Analysis of Principal Components (DAPC) ####
# Transforms data (?), performs PCA, performs Linear Discriminant Analysis (LDA) with PCs
dapc1 <- dapc(my.data, n.pca = 10, n.da = 1) # n.pca = number axes to be retained; n.da = number of axes retained in Discriminant Analysis step
dapc1 # Provides information such as proportion conserved variance and object details

# Density plot of samples along discriminant function 1
pdf(file = paste0(output.dir, "dapc_all_pops.pdf"), width = 10.5, height = 7)
scatter(dapc1, scree.da = F, bg = "white", legend = T
        , txt.leg=rownames(dapc1$means)
        , posi.leg = "topright"
        #, col = my.cols
        )
dev.off()

# assignplot hasn't been explored/implemented here yet
# plot showing the probababilities of assignment of individuals to the different clusters
# assignplot(x = dapc1)

## compoplot not currently working
# ## Composition plot (barplot showing the probabilities of assignments of individuals to the different clusters)
# # compoplot requires a specific library
# install_version("hierfstat", version = "0.04-22", repos = "http://cran.us.r-project.org") # compoplot functional
# library(hierfstat)
# 
# install_version("adegenet", version = "2.0.1", repos = "http://cran.us.r-project.org") # compoplot functional
# library("adegenet")
# 
# par(mar=c(10,3,3,3))
# pdf(file = paste0(output.dir, "compoplot_all_samples.pdf"), width = 24, height = 8)
# compoplot(dapc1
#           #, lab="" # comment out if you want sample labels
#           , txt.leg = rownames(dapc1$means)
#           , posi = "topright"
# #          , cex = 0.7
#           , cex.names = 0.6 
#           )
# dev.off()

# Loading plot # Plot marker variance contribution to DAPC
par(mfrow=c(1,1), mar=c(3,4,3,3))

# Plot the loading values of the different markers into the DAPC
loadingplot(dapc1$var.contr, thres=1e-3)

# # Save out the loading values for the dapc linked to the marker names
# par(mar=c(5,5,4,4))
# dapc1$pca.loadings[1:5,]
# max(dapc1$pca.loadings[,1])
# max(dapc1$var.contr) # this is what's being plotted
# plot(dapc1$pca.loadings[,1], dapc1$var.contr ) # absolute value of pca loadings are 
# plot(abs(dapc1$pca.loadings[,1]), dapc1$var.contr ) # absolute value of pca loadings are highly correlated to the 'var.contrib'
# # Where are the locus names? They aren't included in this, which seems to just take the 'index' of the locus. 
# # Therefore, we have to currently assume that the order is the same as in my.data
# locus.and.dapc1.var.contrib <- as.data.frame(cbind(locNames(my.data), dapc1$var.contr))
# colnames(locus.and.dapc1.var.contrib) <- c("mname","dapc1.var.contr")
# write.csv(x = locus.and.dapc1.var.contrib, file = "11-other_stats/locus_and_dapc1_var_contrib.csv"
#           , quote = F, row.names = F)


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

# Generate population attributes for genind
pop(my.data.gid) # currently null
pop(my.data.gid) <- pop(my.data) # use pop from the my.data object
# or alternately:
# pop(my.data.gid) <- gsub(x = indNames(my.data.gid), pattern = "\\_.*", replacement = "") # create pop ID

# Show sample size per population
table(pop(my.data.gid))
pdf(file = paste0(output.dir, "sample_size.pdf"), width = 8, height = 5)
par(mfrow=c(1,1), mar=c(8,5,3,3))
barplot(table(pop(my.data.gid)), col=funky(17)
        #, las=3, las = 1
        , las=2
        , xlab=""
        , ylab="Sample size"
        , ylim = c(0,40))
#abline(h = c(10,20,30), lty=2)
dev.off()

# Change from genind to hierfstat
all.data.hf <- genind2hierfstat(my.data.gid)
rownames(all.data.hf) <- indNames(my.data.gid)


##### Pairwise Fst #####
pairwise.wc.fst <- pairwise.WCfst(all.data.hf)
write.csv(x = pairwise.wc.fst, file = paste0(output.dir, "all_data_wcfst.csv"))

# Bootstrapping
nboots <- 1000
# requires latest hierfstat (v0.04-29) otherwise get error
# library(devtools)
# install_github("jgx65/hierfstat")
# library("hierfstat")

# boot.fst.all <- boot.ppfst(dat = all.data.hf, nboot = 1000, quant = c(0.025,0.975))
# boot.fst.all
# # note that nboot = 1000 is about the same as nboot = 10,000 (very marginally different)
# 
# # Collect output
# lower.limit <- t(boot.fst.all$ll)
# upper.limit <- boot.fst.all$ul
# upper.limit[is.na(upper.limit)] <- 0
# lower.limit[is.na(lower.limit)] <- 0
# boot.fst.all.output <- upper.limit + lower.limit
# boot.fst.all.output
# 
# filename <- paste0("11-other_stats/all_data_fst_nboot_", nboots, ".csv")
# write.csv(x = boot.fst.all.output, file = filename)

# Heatmap output
# install.packages("plsgenomics")
# require(plsgenomics)
# matrix.heatmap(mat = boot.fst.all.output)

## Good heatmaps here
# require("lattice")
# levelplot(boot.fst.all.output, scales=list(x=list(rot=90)))
# 
# require("gplots")
# heatmap.2(boot.fst.all.output)


save.image(file = paste0(output.dir, "adegenet_output.RData")) # Save out this data

# #### Other unused steps ####
# # Descriptive Statistics on genind file
# div <- summary(my.data.gid) # genetic diversity w/ adegenet
# str(div) # What is contained within this summary object
# # contains: sample size, sample size by pop, number loci and names, population IDs, % NA, Hobs, Hexp
# 
# # Provide three F statistics (Fst (pop/total), Fit (Ind/total), Fis (ind/pop))
# fstat(my.data.gid)
# 
# # which markers are greater than a specific Fst value?
# val <- 0.1
# rownames(basicstat$perloc)[which(basicstat$perloc[, "Fst"] > val)]
# high.fst.markers <- rownames(basicstat$perloc)[which(basicstat$perloc[, "Fst"] > val)]
# length(high.fst.markers)
# 
# # extract these from the genlight object
# head(locNames(my.data.gid))
# high.fst.markers.true <- gsub(x = high.fst.markers, pattern = "^X", replacement = "", perl = T)
# 
# my.data.high.fst <- my.data.gid[, loc=high.fst.markers.true]
# my.data.high.fst
# 
# 
