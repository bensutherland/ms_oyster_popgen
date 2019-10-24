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


# #### INVESTIGATE SCOREPLOT ### WORKING HERE ####
# # get sample names
# file_of_interest <- files[2]
# file_of_interest
# 
# my_vcf <- read.vcfR(file = file_of_interest)
# sample_names <- colnames(my_vcf@gt)[2:length(colnames(my_vcf@gt))]
# pop <- gsub(pattern = "_.*", replacement = "", x = sample_names)
# 
# pcadapt::score_plot(x = pcadapt.out.list[[2]], i = 1, j = 2, pop = pop) #DPB_PEN
# pcadapt::score_plot(x = pcadapt.out.list[[2]], i = 3, j = 4, pop = pop) #DPB_PEN
# pcadapt::score_plot(x = pcadapt.out.list[[2]], i = 5, j = 6, pop = pop) #DPB_PEN
# pcadapt::score_plot(x = pcadapt.out.list[[2]], i = 7, j = 8, pop = pop) #DPB_PEN
# 
# #### END WORKING HERE ####


##### Plot screeplot and scoreplot to ID the number and identity of PCs to obtain #####
contrast.of.interest <- NULL; selected_K <- NULL ; selected_K.list <- list()
my_vcf.FN <- NULL; my_vcf <- NULL; sample_names <- NULL; pop <- NULL

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



##### Manual inspection of results and input of k and PCs of relevance #####

# Now inspect the screeplot to find out how many PCs actually explain a significant proportion of the variation 
# (i.e. to the left of the plateau)

### Manually input the k selection based on the screeplots produced above:
names(pcadapt.out.list)
selected_K.list[["QDC_CHN"]] <- 7
selected_K.list[["GUR_FRA"]] <- 9
selected_K.list[["RSC_CHN"]] <- 4
selected_K.list[["ROS_PIP"]] <- 5
selected_K.list[["DPB_PEN"]] <- 4

selected_K.list[["PEN_PENF"]] <- 1
selected_K.list[["FRA_FRAF"]] <- 1
selected_K.list[["CHN_CHNF"]] <- 1


# Now inspect the score plots to find out which of the PCs within the k selection above 
# actually separate the groups of interest
# selected_PCs.list[["QDC_CHN"]] <- 
# selected_PCs.list[["GUR_FRA"]] <- 
# selected_PCs.list[["RSC_CHN"]] <- 
# selected_PCs.list[["ROS_PIP"]] <- 
# selected_PCs.list[["DPB_PEN"]] <- 
# 
# selected_PCs.list[["PEN_PENF"]] <- 
# selected_PCs.list[["FRA_FRAF"]] <- 
# selected_PCs.list[["CHN_CHNF"]] <- 




#### Run pcadapt for each contrast including the PCs as selected above manually
# Using the input obj pcadapt.obj.list and selected_K.list, re-run pcadapt using dynamically selected k
pcadapt.out_dynamic_K.list <- list()

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


#### Choosing cutoff for outlier detection ####
qvals <- NULL; outliers.list <- list(); pval_and_qval.df <- NULL; num_outliers.list <- list(); mname.df <- NULL

# Extract qvals
for(i in 1:length(pcadapt.out_dynamic_K.list)){
  
  # Name for round
  contrast.of.interest <- names(pcadapt.out_dynamic_K.list)[i]
  
  # extract pvals and convert to qval
  qvals <- qvalue(pcadapt.out_dynamic_K.list[[contrast.of.interest]]$pvalues)$qvalues
  outliers.list[[contrast.of.interest]] <- qvals
  
  # Combine pval with qvalue
  pval_and_qval.df <- as.data.frame(cbind(pcadapt.out_dynamic_K.list[[contrast.of.interest]]$pvalues, qvals), stringsAsFactors = F)
  colnames(pval_and_qval.df) <- c("pval", "qval")
  
  # Bring in marker name
  mname.df <- read.table(file = paste0("13_selection/", contrast.of.interest, "_mnames.txt"), header = F, sep = "\t")
  mname.df$V2 <- seq(1:length(mname.df$V1))
  colnames(mname.df) <- c("mname", "index")
  
  # TODO: add test to make sure the two df match before combining them
  
  # Combine pval/qval with marker name
  pval_and_qval.df <- cbind(mname.df, pval_and_qval.df)
  
  write.table(file = paste0("13_selection/", "qvals_", contrast.of.interest, ".csv")
              , x = pval_and_qval.df, sep = ",", col.names = T, row.names = F
              , quote = F)
  
  # # Finds the *indices* (not names) of which markers are outliers
  num_outliers.list[[contrast.of.interest]] <- length(which(qvals < qval_cutoff))
  
}

num_outliers.list

num_outliers <- sapply(num_outliers.list, function(x){as.numeric(x[1])})
as.data.frame(num_outliers)
num_outliers.df <- as.data.frame(as.numeric(num_outliers), stringsAsFactors = FALSE)
num_outliers.df$contrast <- as.character(names(num_outliers))
colnames(num_outliers.df) <- c("outliers", "contrast")
num_outliers.df <- num_outliers.df[, c("contrast", "outliers")]

write.csv(x = num_outliers.df, file = paste0("13_selection/number_outliers_per_contrast_qval_", qval_cutoff, ".csv"), quote = F, row.names = F)
