# Separate Analysis of individual populations
# Requires the running of 'adegenet.R'

#### Front Matter ####
# Clean space
# rm(list=ls())

# Load required libraries
library("purrr")
library("dplyr")
library("ape")
library("pegas")
library("seqinr")
library("ggplot2")
library("devtools")
library("hierfstat")
library("adegenet")
library("lattice")
library("gplots")

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

#### 01. Input data and prepare ####
## Load inputs
# Load part 1 results
load(file = paste0(output.dir, "adegenet_output.RData"))
my.data.gid

# Load colours file
my_cols.df <- read.csv(file = "../ms_oyster_popgen/00_archive/my_cols.csv", stringsAsFactors = F)
str(my_cols.df)

# Load sample phenotypes file (i.e. for Size_Class)
sample.data <- read.csv(file = "../ms_oyster_popgen/00_archive/master_list_all_samples.csv", header = T)
head(sample.data)
# Keep only the necessary columns
sample.data <- sample.data[,c("Sample_Number", "Date_Collected", "CollectionSite", "Collection_Site_Location_Details", "Year_Class", "Size_Class")]
head(sample.data)

## Prepare individual analyses of genind by separating all pops
sep.obj <- seppop(x = my.data.gid)
names(sep.obj)

# # The following is to input markers that are not in HWP
# # Load markers not tested for HWE for selection analysis
# load("24-selection/my_sel_data.Rdata")
# my.sel.data.gid
# sep.sel.obj <- seppop(x = my.sel.data.gid)
# # End section

# Create list of various combinations of repooled pops
datatype.list <- list()
datatype.list[["global"]] <- repool(sep.obj$PEN
                                     , sep.obj$CHN, sep.obj$QDC, sep.obj$RSC
                                     , sep.obj$DPB, sep.obj$ROS, sep.obj$GUR
                                     , sep.obj$FRA
                                     , sep.obj$JPN
                                     )
datatype.list[["bc_all"]] <- repool(sep.obj$HIS, sep.obj$PEN, sep.obj$PIP, sep.obj$SER)

datatype.list[["bc_size"]] <- repool(sep.obj$PEN, sep.obj$HIS)

datatype.list[["ch_all"]] <- repool(sep.obj$CHN, sep.obj$CHNF, sep.obj$QDC, sep.obj$RSC)

# and selection items (w/ loci out of HWE still included):
datatype.list[["bc_sel"]] <- repool(sep.obj$PEN, sep.obj$PENF)

datatype.list[["fr_sel"]] <- repool(sep.obj$FRA, sep.obj$FRAF)

datatype.list[["ch_sel"]] <- repool(sep.obj$CHN, sep.obj$CHNF)

datatype.list[["dpb_sel"]] <- repool(sep.obj$DPB, sep.obj$PEN)

datatype.list[["ros_sel"]] <- repool(sep.obj$ROS, sep.obj$PIP)

datatype.list[["qdc_sel"]] <- repool(sep.obj$QDC, sep.obj$CHN)

datatype.list[["rsc_sel"]] <- repool(sep.obj$RSC, sep.obj$CHN)

datatype.list[["gur_sel"]] <- repool(sep.obj$GUR, sep.obj$FRA)

# rm(sep.obj)


#### 02. Analyze each datatype ####

# Loop over all datatypes
pops.involved <- NULL; subset.cols <- NULL

for(i in 1:length(datatype.list)){
  
  # Set datatype and select from the list
  datatype <- names(datatype.list)[i]
  data.gid <- datatype.list[[datatype]]
  print(paste("Now working on datatype", datatype))
  
  # Get the colors for this datatype
  pops.involved <- levels(datatype.list[[datatype]]$pop) # identify pops in this datatype
  pops.involved.df <- as.data.frame(pops.involved) # make the vector of pops involved a df so that can use merge
  my_cols_per_datatype.df <- merge(x = pops.involved.df, y = my_cols.df, by.x = "pops.involved", by.y = "my.pops", sort = F)
  cols <- my_cols_per_datatype.df[,"my.cols"] # note: order is based on the order of levels
  
  # Set your data into the data.gid object
  data.gid <- datatype.list[[datatype]]
  nInd(data.gid)  

  # Incorporate size phenotype (used only for bc_size)
  # Obtain sample numbers to match w/ sample.info database
  indiv.used <- gsub(pattern = "*.*\\_", replacement = "", indNames(data.gid))
  # indiv.used <- gsub(pattern = "\\..*", replacement = "", indiv.used.pre) # Can remove now?
  indiv.used <- as.numeric(indiv.used)

  # Subset sample.data for the samples in this datatype
  subset.of.sample.data <- sample.data[indiv.used, c("Sample_Number", "CollectionSite", "Size_Class")]
  head(subset.of.sample.data)

  # Only keep the first three letters of the location for short name
  subset.of.sample.data$CollectionSite <- substr(subset.of.sample.data$CollectionSite, start = 1, stop = 3)
  subset.of.sample.data$CollectionSite <- toupper(subset.of.sample.data$CollectionSite)
  head(subset.of.sample.data)

  # Combine the CollectionSize and Size_Class cols together
  subset.of.sample.data$pop.id <- apply(subset.of.sample.data[, c("CollectionSite","Size_Class")], 1, paste, collapse = "_")
  subset.of.sample.data$pop.id
  head(subset.of.sample.data)

  # Confirm lengths match
  if(length(indNames(data.gid)) == length(subset.of.sample.data$pop.id)){ 
    print("All samples appear to be accounted for")
  } else {
      print("There are some samples that haven't been assigned pop IDs, missing from db?")}

  # For the datatype that uses a size analysis, replace pop character with the one that includes size info
  if(datatype == "bc_size"){ 
    print("**Replacing population ID to include size data**")
    pop(data.gid) <- subset.of.sample.data$pop.id
    
    # fix cols
    cols <- c(rep(x = cols[1], times = 2), rep(x = cols[2], times = 2))
    
  } else {
    print("No size class information required")
  }

  pop(data.gid)

  # Identify the populations in the data
  levels(pop(x = data.gid))
  nPop(data.gid)
  nInd(data.gid)
  indNames(data.gid)


  # Show sample size per population
  filename <- paste(output.dir, datatype, "_sample_size_per_pop.pdf", sep = "")
  pdf(file = filename, width = 7, height = 4)
  par(mfrow=c(1,1), mar=c(8,5,3,3))
  barplot(table(pop(data.gid))
          , col=cols
        , las=2
        , xlab=""
        , ylab="Sample size"
        , ylim = c(0,40)
        )
  abline(h = c(10,20,30), lty=2)
  dev.off()

  # Change from genind file to hfstat
  data.hf <- genind2hierfstat(data.gid)
  rownames(data.hf) <- indNames(data.gid)

  ## PCA ## on a matrix of individual genotype frequencies (hierfstat)
  y <- indpca(data.hf, ind.labels = rownames(data.hf))
  #y <- indpca(data.hf, ind.labels = pop(data.gid)) # this allows to view the size class type, if wanted

  filename <- paste(output.dir, datatype, "_sample_PCA.pdf", sep = "")
  pdf(file = filename, width = 11, height = 6)
  plot(y, cex = 0.5)
  dev.off()

  ## DAPC ##
  dapc <- dapc(data.gid, n.pca = 10, n.da = 1)

  filename <- paste(output.dir, datatype, "_sample_DAPC.pdf", sep = "")
  pdf(file = filename, width = 11, height = 6)
  scatter(dapc, scree.da = F, bg = "white", legend = T
        , txt.leg=rownames(dapc$means)
        , posi.leg = "topleft"
        , col = cols
        )
  dev.off()

  # Loading plot to plot marker variance contribution to DAPC
  filename <- paste(output.dir, datatype, "_DAPC_loadings.pdf", sep = "")
  pdf(file = filename, width = 11, height = 4)
  par(mfrow=c(1,1), mar=c(3,4,3,3))
  loadingplot(dapc$var.contr
              , thres=2e-3
              , las = 1
              #, lab.jitter =1
              #, srt = 90 # makes text vertical
              )
  dev.off()

  ## FST Analysis ##
  # Pairwise Fst (hierfstat)
  pairwise.wc.fst <- pairwise.WCfst(data.hf)

  filename <- paste(output.dir, datatype, "_pairwise_wc_fst.csv", sep = "")
  write.csv(pairwise.wc.fst, file = filename)

  # Pairwise Fst w/ bootstrapping (hierfstat)
  boot.fst <- boot.ppfst(dat = data.hf, nboot = 1000, quant = c(0.025,0.975))
  boot.fst

  # Collect output
  lower.limit <- t(boot.fst$ll)
  upper.limit <- boot.fst$ul
  upper.limit[is.na(upper.limit)] <- 0
  lower.limit[is.na(lower.limit)] <- 0
  boot.fst.output <- upper.limit + lower.limit
  boot.fst.output

  filename <- paste(output.dir, datatype, "_boot_fst_output.csv", sep = "")
  write.csv(x = boot.fst.output, file = filename)

  # Plot Fst results in heatmap
  filename <- paste(output.dir, datatype, "_levelplot_fst_heatmap.pdf", sep = "")
  pdf(file = filename, width = 10.5, height = 6)
  
  levelplot(boot.fst.output, scales=list(x=list(rot=90))
          , xlab = "Fst (lower limit Fst)"
          , ylab = "Fst (upper limit Fst)"
          )
  dev.off()

  filename <- paste(output.dir, datatype, "_heatmap2_fst_heatmap.pdf", sep = "")
  pdf(file = filename, width = 6, height = 6)
  heatmap.2(boot.fst.output)
  dev.off()

  #### Per loc stats ####
  # Basic stats (gives per loc and overall, but only one value, not one for each comparison)
  perloc.stats <- basic.stats(data.hf, diploid = T, digits = 4)
  str(perloc.stats$perloc) # this is the key per locus data from this analysis

  # Per loc Fst figure
  filename <- paste(output.dir, datatype, "_per_loc_overall_fst.pdf", sep = "")
  pdf(file = filename, width = 9, height = 6)
  par(mfrow=c(1,1), mar=c(5,5,2.5,3))
  plot(perloc.stats$perloc$Fst, ylab = "Fst", las = 1) # Plot Fst
  text(x = 2000, y = 0.2, labels = paste(ncol(data.hf)-1,"loci" ))
  text(x = 2000, y = 0.19, labels = paste(nrow(data.hf),"ind" ))

  dev.off()


  # Save out these Fst values per locus to match to location in the genome
  perloc.stats.df <- perloc.stats$perloc # save df w/ Ho, Hs, Ht, Dst, Htp, Dstp, Fst, Fstp, Fis, Dest (from hierfstat)
  filename <- paste(output.dir, datatype, "_perloc_stats.csv", sep = "")
  write.csv(x = perloc.stats.df, file = filename)

}
