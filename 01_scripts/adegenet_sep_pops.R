# Separate Analysis of individual populations

#### Front Matter ####
# Clean space
# rm(list=ls())

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

# Set working directory
# Xavier
setwd("~/Documents/01_moore_oyster_project/stacks_workflow_no_Qc")

# Load part 1 results
load("11-other_stats/adegenet_output.RData")

my.data.gid

# Read in database with sample phenotypes, most notably the Size_Class attribute
sample.data <- read.csv(file = "../ms_oyster_popgen/00_archive/master_list_all_samples.csv", header = T)
head(sample.data)
# reduce
sample.data <- sample.data[,c("Sample_Number", "Date_Collected", "CollectionSite", "Collection_Site_Location_Details", "Year_Class", "Size_Class")]
head(sample.data)


# Separate gid
sep.obj <- seppop(x = my.data.gid)
names(sep.obj)

# Create a subsetted list of repooled subcomponents of the full data
datatype.list <- list()
datatype.list[["all.gid"]] <- repool(sep.obj$Hisnit, sep.obj$Pendrell, sep.obj$PendrellFarm, sep.obj$Pipestem, sep.obj$Serpentine
                                     , sep.obj$ChinaBe, sep.obj$ChinaQDC, sep.obj$ChinaQDW, sep.obj$ChinaRSC
                                     , sep.obj$DeepBay, sep.obj$Rosewall, sep.obj$Guernsey
                                     , sep.obj$FranceW, sep.obj$FranceC
                                     )
datatype.list[["all.bc.gid"]] <- repool(sep.obj$Hisnit, sep.obj$Pendrell, sep.obj$PendrellFarm, sep.obj$Pipestem, sep.obj$Serpentine)

datatype.list[["bc.wild.size.gid"]] <- repool(sep.obj$Serpentine, sep.obj$Hisnit) # requires use of sample.data for phenotype of size

datatype.list[["bc.sel.gid"]] <- repool(sep.obj$Pendrell, sep.obj$PendrellFarm) # PendrellFarm is culture selection
datatype.list[["fr.gid"]] <- repool(sep.obj$FranceW, sep.obj$FranceC) # FranceC is culture selection
datatype.list[["ch.gid"]] <- repool(sep.obj$ChinaBe, sep.obj$ChinaQDW) # QDW is culture selection


###
# Several notes
# ChinaQDW = spat taken from culture region (QD) by artificial attachment, then grown in 'lantern net culture'
# ChinaBe = from QD region, sampled directly from rock on the beach
###

#### Choose dataset ####
datatype <- "all.gid"
# datatype <- "all.bc.gid"

# datatype <- "bc.wild.size.gid"

# datatype <- "bc.sel.gid"
# datatype <- "fr.gid"
# datatype <- "ch.gid"

# Set your data into the data.gid object
data.gid <- datatype.list[[datatype]]
  

#### Analysis ####

# Include size data if needed (note: can run if not needed, will not replace popid
# Obtain sample numbers only (no pop info) to match w/ sample.info
indiv.used <- as.numeric(gsub(pattern = "*.*\\_", replacement = "", indNames(data.gid)))

# Take the samples from the sample.data, specifically with the info on the location and size class
subset.of.sample.data <- sample.data[indiv.used, c("CollectionSite", "Size_Class")]
# Only keep the first three letters of the location
subset.of.sample.data$CollectionSite <- substr(subset.of.sample.data$CollectionSite, start = 1, stop = 3)

# collapse the columns together
subset.of.sample.data$pop.id <- apply(subset.of.sample.data[,c(1:ncol(subset.of.sample.data))], 1, paste, collapse = "_")
subset.of.sample.data$pop.id

# Confirm lengths match
if(length(indNames(data.gid)) == length(subset.of.sample.data$pop.id)){ print("yes")} else {print("no")}

# For only the datatype bc.wild.size.gid, replace the pop id with the one including oyster sizes
if(datatype == "bc.wild.size.gid"){ 
  print("Replacing population ID to include size data")
  pop(data.gid) <- subset.of.sample.data$pop.id
  print(pop(data.gid))
} else {
    print("No size class information required")
  }

pop(data.gid)

# Identify the populations in the data
levels(pop(x = data.gid))
nPop(data.gid)
nInd(data.gid)
indNames(data.gid)
pop(data.gid)

# Show sample size per population
filename <- paste("11-other_stats/", datatype, "_sample_size_per_pop.pdf", sep = "")
pdf(file = filename, width = 7, height = 4)
par(mfrow=c(1,1), mar=c(8,5,3,3))
barplot(table(pop(data.gid)), col=funky(17)
        , las=2
        , xlab=""
        , ylab="Sample size"
        , ylim = c(0,40))
abline(h = c(10,20,30), lty=2)
dev.off()

# Change from genind file to hfstat
data.hf <- genind2hierfstat(data.gid)
rownames(data.hf) <- indNames(data.gid)

# PCA on a matrix of individual genotype frequencies (hierfstat)
y <- indpca(data.hf, ind.labels = rownames(data.hf))
#y <- indpca(data.hf, ind.labels = pop(data.gid)) # this allows to view the size class type, if wanted

filename <- paste("11-other_stats/", datatype, "_sample_PCA.pdf", sep = "")
pdf(file = filename, width = 11, height = 6)
plot(y, cex = 0.7)
dev.off()

# DAPC
dapc <- dapc(data.gid, n.pca = 10, n.da = 1)

filename <- paste("11-other_stats/", datatype, "_sample_DAPC.pdf", sep = "")
pdf(file = filename, width = 11, height = 6)
scatter(dapc, scree.da = F, bg = "white", legend = T
        , txt.leg=rownames(dapc$means)
        , posi.leg = "topleft")
dev.off()

# Composition plot (barplot showing the probabilities of assignments of individuals to the different clusters)
filename <- paste("11-other_stats/", datatype, "_compoplot.pdf", sep = "")
pdf(file = filename, width = 11, height = 6)
par(mar=c(10,3,3,3))
compoplot(dapc
          , txt.leg = rownames(dapc$means)
          , posi = "topright"
          , cex.names = 0.6 
)
dev.off()

# Loading plot # Plot marker variance contribution to DAPC
filename <- paste("11-other_stats/", datatype, "_DAPC_loadings.pdf", sep = "")
pdf(file = filename, width = 11, height = 4)
par(mfrow=c(1,1), mar=c(3,4,3,3))
loadingplot(dapc$var.contr, thres=1e-3, las = 1)
dev.off()

# Pairwise Fst (hierfstat)
pairwise.wc.fst <- pairwise.WCfst(data.hf)

filename <- paste("11-other_stats/", datatype, "_pairwise_wc_fst.csv", sep = "")
write.csv(pairwise.wc.fst, file = filename)

# Pairwise Fst w/ bootstrapping (hierfstat)
# requires latest hierfstat (v0.04-29) otherwise get error
# library(devtools)
# install_github("jgx65/hierfstat")
library("hierfstat")
boot.fst <- boot.ppfst(dat = data.hf, nboot = 1000, quant = c(0.025,0.975))
boot.fst

# Collect output
lower.limit <- t(boot.fst$ll)
upper.limit <- boot.fst$ul
upper.limit[is.na(upper.limit)] <- 0
lower.limit[is.na(lower.limit)] <- 0
boot.fst.output <- upper.limit + lower.limit
boot.fst.output

filename <- paste("11-other_stats/", datatype, "_boot_fst_output.csv", sep = "")
write.csv(x = boot.fst.output, file = filename)

# Heatmap output
# install.packages("plsgenomics")
# require(plsgenomics)
# matrix.heatmap(mat = boot.fst.all.output)

filename <- paste("11-other_stats/", datatype, "_levelplot_fst_heatmap.pdf", sep = "")
pdf(file = filename, width = 10.5, height = 6)
require("lattice")
levelplot(boot.fst.output, scales=list(x=list(rot=90))
          , xlab = "Fst (lower limit Fst)"
          , ylab = "Fst (upper limit Fst)")
dev.off()

require("gplots")
heatmap.2(boot.fst.output)

## Per loc Fst
# Basic stats (gives per loc and overall, but only one value, not one for each comparison)
perloc.stats <- basic.stats(data.hf, diploid = T, digits = 4)
str(perloc.stats$perloc) # this is the key data from this analysis

# Per loc Fst figure
filename <- paste("11-other_stats/", datatype, "_per_loc_overall_fst.pdf", sep = "")
pdf(file = filename, width = 9, height = 6)
par(mfrow=c(1,1), mar=c(5,5,2.5,3))
plot(perloc.stats$perloc$Fst, ylab = "Fst", las = 1) # Plot Fst
text(x = 2000, y = 0.2, labels = paste(ncol(data.hf)-1,"loci" ))
text(x = 2000, y = 0.19, labels = paste(nrow(data.hf),"ind" ))

dev.off()


# Save out these Fst values per locus to match to location in the genome
perloc.stats.df <- perloc.stats$perloc # save df w/ Ho, Hs, Ht, Dst, Htp, Dstp, Fst, Fstp, Fis, Dest (from hierfstat)
filename <- paste("11-other_stats/", datatype, "_perloc_stats.csv", sep = "")
write.csv(x = perloc.stats.df, file = filename)
# match these to their genomic location and plot in Rscript (#todo)



