# Method to detect outliers, developed by Anne-Laure Ferchaud (U Laval) and Ben Sutherland (DFO)
# Uses pcadapt methods as described here: https://cran.r-project.org/web/packages/pcadapt/vignettes/pcadapt.html

#### Front Matter ####
# Clean space
# rm(list=ls())

# Set working directory for stark or xavier
# if on a different system, prompt for working directory
if(Sys.info()["nodename"] == "stark"){ 
  print("On Stark, ready to go")
  setwd("/mnt/data/01_moore_oyster_project/stacks_workflow/") # stark
} else if(Sys.info()["nodename"] == "Xavier"){
  print("On Xavier, ready to go")
  setwd("/hdd/01_moore_oyster_project/stacks_workflow/") # Xavier
} else {
  print("You are on an unrecognized system, please set working directory manually")
}

## Info
# sessionInfo()

### Load libraries
# install.packages("pcadapt")
# source("https://bioconductor.org/biocLite.R")
# biocLite("qvalue")

library(pcadapt)
library(qvalue)

# Set variables
qval_cutoff <- 0.01

# Identify files
files <- list.files(path = "13_selection/", pattern = "recode.vcf", full.names = T)

#### Run pcadapt for each contrast with full pcs retained ####
# This is to determine which PCs contribute significantly to the data, and 
# to determine which PCs separate the groups
# set nulls
pcadapt.obj <- NULL; pcadapt.FN <- NULL; pcadapt.out <- NULL
pcadapt.out.list <- list(); obj.name <- NULL; pcadapt.obj.list <- list()

for(i in 1:length(files)){
  
  # Select file
  pcadapt.FN <- files[i]
  print(paste0("Working on ", pcadapt.FN))
  
  # Create object name from filename (e.g. ROS_PIP)
  obj.name <- gsub(pattern = ".*//", replacement = "", x = files[i])
  obj.name <- gsub(pattern = "\\..*", replacement = "", x = obj.name)
  
  # Read in data
  pcadapt.obj <- read.pcadapt(input = pcadapt.FN, type = "vcf")
  
  # Save input data into a list using the object name above
  pcadapt.obj.list[[obj.name]] <- pcadapt.obj
  
  # Run pcadapt
  print("Running pcadapt with K = 20")
  pcadapt.out <- pcadapt(input = pcadapt.obj, K = 20)
  
  # Save output
  pcadapt.out.list[[obj.name]] <- pcadapt.out
  
}

# Note: the pcadapt.obj.list contains all the input data


# Example of what output data looks like:
# View the objects in the pcadapt list
summary(pcadapt.out.list[[1]])

# Use vector containing the K-ordered, squared root of the proportion of variance explained by each PC
pcadapt.out.list[[1]]$singular.values

# Calculate the proportion of variance explained by each PC (by squaring)
pcadapt.out.list[[1]]$singular.values ^ 2
# End example

##### Plot screeplot and scoreplot to ID the number and identity of PCs to obtain #####
contrast.of.interest <- NULL; selected_K <- NULL ; selected_K.list <- list()
my_vcf.FN <- NULL; my_vcf <- NULL; sample_names <- NULL; pop <- NULL
mnames <- NULL; mnames.list <- list()

for(i in 1:length(pcadapt.out.list)){
  
  # Name for round
  contrast.of.interest <- names(pcadapt.out.list)[i]
  print(contrast.of.interest) # reporting
  
  # SCORE PLOT SECTION
  # Obtain the individual names for this round
  my_vcf.FN <- paste0("13_selection//", contrast.of.interest, ".recode.vcf")
  my_vcf <- read.vcfR(file = my_vcf.FN)
  sample_names <- colnames(my_vcf@gt)[2:length(colnames(my_vcf@gt))] # obtain sample names, which start at column 2
  pop <- gsub(pattern = "_.*", replacement = "", x = sample_names) # obtain pop names, by dropping everything after the underscore
  
  # Obtain marker names for this round
  head(my_vcf@fix)
  
  mnames <- paste0(">CLocus_"
                   , gsub(x = my_vcf@fix[,"ID"], pattern = ":.*", replacement = "")
                   , "-"
                   , my_vcf@fix[,"CHROM"]
                )
  
  mnames.list[[contrast.of.interest]] <- mnames
  
  # Use scoreplot, saving 10 pcs (five plots, 2 PCs each) for evaluating
  for(s in seq(from = 1, to = 10, by = 2)){
    
    # Save out five score plots 
    pdf(file = paste0("13_selection/", contrast.of.interest, "_scoreplot_", s, "_", s+1, ".pdf"), width = 5, height = 5)
    print(
      pcadapt::score_plot(x = pcadapt.out.list[[i]], i = s, j = s + 1, pop = pop)
    )
    dev.off()
  }
  # END SCORE PLOT SECTION
  
  # SCREE PLOT SECTION
  pdf(file = paste0("13_selection/", contrast.of.interest, "_screeplot_K20.pdf"), width = 5, height = 5)
  print(
    plot(pcadapt.out.list[[i]], option = "screeplot", K = 20)
  )
  dev.off()
  
}



##### Manually inspect results and input k and PCs of relevance #####

# Inspect the screeplots to find out how many PCs actually explain a significant proportion of the variation 
# (i.e. to the left of the plateau)

### Manually input the k selection based on the screeplots produced above:
names(pcadapt.out.list)
selected_K.list[["QDC_CHN"]] <- 6
selected_K.list[["GUR_FRA"]] <- 6
selected_K.list[["RSC_CHN"]] <- 3
selected_K.list[["ROS_PIP"]] <- 5
selected_K.list[["DPB_PEN"]] <- 3

selected_K.list[["PEN_PENF"]] <- 1
selected_K.list[["FRA_FRAF"]] <- 1
selected_K.list[["CHN_CHNF"]] <- 1


# Inspect the score plots to find out which of the PCs within the k selection above 
# actually separate the groups of interest
selected_PCs.list <- list()
selected_PCs.list[["QDC_CHN"]] <- c(1,2,3,4)
selected_PCs.list[["GUR_FRA"]] <- c(1)
selected_PCs.list[["RSC_CHN"]] <- c(1)
selected_PCs.list[["ROS_PIP"]] <- c(1,2,3,4)
selected_PCs.list[["DPB_PEN"]] <- c(1,3)

selected_PCs.list[["PEN_PENF"]] <- c(0)
selected_PCs.list[["FRA_FRAF"]] <- c(0)
selected_PCs.list[["CHN_CHNF"]] <- c(0)


#### Run pcadapt for each contrast including the PCs as selected above manually
# Using the input obj pcadapt.obj.list and selected_K.list, re-run pcadapt using dynamically selected k
pcadapt.out_dynamic_K.list <- list() # This will be the OUTPUT required for downstream

for(i in 1:length(pcadapt.obj.list)){
  
  # Name for round
  contrast.of.interest <- names(pcadapt.obj.list)[i]
  
  pcadapt.out_dynamic_K.list[[contrast.of.interest]] <- pcadapt(input = pcadapt.obj.list[[contrast.of.interest]]
                                                                , K = selected_K.list[[contrast.of.interest]]
                                                                , method = "mahalanobis"
                                                                , min.maf = 0.05
                                                                , ploidy = 2)
  
  # Save outputs
  # Manhattan plot
  pdf(file = paste0("13_selection/", contrast.of.interest, "_manhattan.pdf"), width = 5, height = 4)
  print(
      pcadapt::manhattan_plot(x = pcadapt.out_dynamic_K.list[[contrast.of.interest]]
                          , chr.info = NULL # Note: could supply chr info here
                          , snp.info = NULL
                          , plt.pkg = "ggplot"
                          )
  )
  dev.off()
  
  # qqplot
  pdf(file = paste0("13_selection/", contrast.of.interest, "_qqplot.pdf"), width = 5, height = 4)
  print(
    pcadapt::qq_plot(x = pcadapt.out_dynamic_K.list[[contrast.of.interest]])
  )
  dev.off()
  
  # pval hist
  pdf(file = paste0("13_selection/", contrast.of.interest, "_pval_hist.pdf"), width = 5, height = 4)
  print(
    hist(pcadapt.out_dynamic_K.list[[contrast.of.interest]]$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
  )
  dev.off()
  
  # stat distrib
  pdf(file = paste0("13_selection/", contrast.of.interest, "_stat_distrib.pdf"), width = 5, height = 4)
  print(
    plot(pcadapt.out_dynamic_K.list[[contrast.of.interest]], option = "stat.distribution")
  )
  dev.off()
  
}
# The output is now in pcadapt.out_dynamic_K.list
# Plots saved in 13_selection


##### HERE TODAY (Self note) #####
# To edit: use get.pc() on each pcadapt object to find which is the top-loading PC for each marker, and retain only those in the
# retained PCs list (using %in% or similar)
# Also: fix the mname issue somewhere above, probably when importing the VCF file
##### END HERE #####


#### Choosing cutoff for outlier detection ####
qvals <- NULL; outliers.list <- list(); outliers_details.df <- NULL; num_outliers.list <- list()
mname.df <- NULL

# Extract qvals from the pcadapt obj
for(i in 1:length(pcadapt.out_dynamic_K.list)){
  
  # Name for round
  contrast.of.interest <- names(pcadapt.out_dynamic_K.list)[i]
  
  # Extract pvals, convert to qval, create df
  outliers_details.df <- pcadapt.out_dynamic_K.list[[contrast.of.interest]]$pvalues
  outliers_details.df <- as.data.frame(outliers_details.df, stringsAsFactors=FALSE)
  colnames(outliers_details.df) <- "pval"
  head(outliers_details.df)
  outliers_details.df$qval <- qvalue(pcadapt.out_dynamic_K.list[[contrast.of.interest]]$pvalues)$qvalues
  str(outliers_details.df)
  
  # Find out which pc contributed most to make this marker significant
  top.pc <- get.pc(x = pcadapt.out_dynamic_K.list[[contrast.of.interest]], list = seq(1:length(qvals)))
  colnames(top.pc) <- c("index", "top.pc")
  
  # Bring this into the df too
  outliers_details.df <- cbind(outliers_details.df, top.pc)
  head(outliers_details.df)
  
  # And finally add the SNP name generated earlier
  mname.df <- mnames.list[[contrast.of.interest]]
  mname.df <- as.data.frame(mname.df)
  mname.df[,"index"] <- seq(1:length(mname.df$mname.df))
  colnames(mname.df) <- c("mname", "index2")
  head(mname.df)
  outliers_details.df <- cbind(outliers_details.df, mname.df)
  
  head(outliers_details.df)
  
  # view results briefly
  table(outliers_details.df$top.pc)
  table(outliers_details.df$top.pc[outliers_details.df$qval < qval_cutoff])
  
  # Write out results
  write.table(file = paste0("13_selection/", "qvals_", contrast.of.interest, ".csv")
              , x = outliers_details.df, sep = ",", col.names = T, row.names = F
              , quote = F)
  
  # Find the number of markers that are outliers for the pcs of interest
  num_outliers.list[[contrast.of.interest]] <- table(outliers_details.df$top.pc %in% selected_PCs.list[[contrast.of.interest]])["TRUE"]
}

num_outliers.list

# 
num_outliers <- sapply(num_outliers.list, function(x){as.numeric(x[1])})
as.data.frame(num_outliers)
num_outliers.df <- as.data.frame(as.numeric(num_outliers), stringsAsFactors = FALSE)
num_outliers.df$contrast <- as.character(names(num_outliers))
colnames(num_outliers.df) <- c("outliers", "contrast")
num_outliers.df <- num_outliers.df[, c("contrast", "outliers")]

write.csv(x = num_outliers.df, file = paste0("13_selection/number_outliers_per_contrast_qval_", qval_cutoff, ".csv"), quote = F, row.names = F)
