# Method to detect outliers, developed by Anne-Laure Ferchaud (U Laval) and Ben Sutherland (DFO)
# Uses RDA methods as per <to update>

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

## Install libraries 
# Note left in the code: replace "-1" by "NA"
library(psych)
library(vegan)

# Set datatypes
datatypes <- c("FARM", "DOMESTICATION")


#### 01. Read in data ####

for(i in 1:length(datatypes)){
  
  # Select your datatype
  datatype <- datatypes[i]
  print(paste0("Working on ", datatype))
  
  
  ##  Read in genotype data ##
  geno.FN <- paste0("13_selection/rda_analysis/", datatype, "_geno.012")
  ind.FN <- paste0("13_selection/rda_analysis/", datatype, "_geno.012.indv")
  snp.FN <- paste0("13_selection/rda_analysis/", datatype, "_geno.012.pos")
  
  geno <- read.table(file = geno.FN, sep = "\t")
  geno[1:5,1:5]
  geno <- geno[, -1] # drop the index column (counts from 0)
  
  ## Read in the individuals data and use as rownames for the geno obj
  ind <- read.table(file = ind.FN)
  row.names(geno) <- ind$V1
  geno[1:5,1:5]

  ## Read in the marker names and use as column names
  snp <- read.table(file = snp.FN)
  snp <- paste(snp$V1, snp$V2, sep = "_")
  colnames(geno) <- snp
  geno[1:5,1:5]
  
  # Convert to matrix
  gen <- as.matrix(geno)
  gen[1:5,1:5]
  
  
  
  
  # Identify NAs
  sum(is.na(gen)) # no NAs
  gen[gen == -1] <- NA # when using --012 format from VCFtools, NA are designated -1
  sum(is.na(gen)) # NAs now in data
  
  # Impute
  gen.imp <- apply(X = gen, MARGIN = 2, function(x)replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
  # work over columns (replace the missing data with the most common genotype in the dataset (?))
  sum(is.na(gen.imp))
  
  #### 02. Make pop level data (env.df) ####

  # e.g. in the following format:
  # INDIVIDUALS	POP	FARM
  # CHN_1551	CHN	WILD
  # CHN_1552	CHN	WILD

  ### OLD METHOD ###
  # Read in the samples in this comparison
  # env.FN <- paste0("13_selection/rda_analysis/", datatype, ".txt")
  # env <- read.table(env.FN)
  ### END OLD METHOD ###
  
  # Use the gen object itself to identify individuals
  INDIVIDUALS <- rownames(gen.imp)
  env.df <- as.data.frame(INDIVIDUALS)
  env.df$POP <- gsub(pattern = "_.*", replacement = "", x =  env.df$INDIVIDUALS)
  head(env.df)
  
  # Identify which pops are farm in either dataset
  if(datatype=="FARM"){
    # Identify the farm pops
    selection_pops <- c("CHNF", "PENF", "FRAF")
  }else if(datatype=="DOMESTICATION"){
    # Identify the farm pops
    selection_pops <- c("DPB", "GUR", "QDC", "ROS", "RSC")
  }
  
  # Identify which are the selection pops
  env.df$FARM[env.df$POP %in% selection_pops] <- "FARM"
  env.df$FARM[!env.df$POP %in% selection_pops] <- "WILD" # find opposite of the above
  head(env.df)
  tail(env.df)
  
  # Bring in colour info
  my_cols.df <- read.csv(file = "../ms_oyster_popgen/00_archive/my_cols.csv", stringsAsFactors = F) 
  str(my_cols.df)
  
  # Combine colour info
  env.df <- merge(x = env.df, y = my_cols.df, by.x = "POP", by.y = "my.pops", sort = F, all.x = T)
  head(env.df)
  
  # Create shape to plot
  env.df$shape <- env.df$FARM
  env.df$shape <- gsub(pattern = "FARM", replacement = "18", x = env.df$shape)
  env.df$shape <- gsub(pattern = "WILD", replacement = "17", x = env.df$shape)
  env.df$shape <- as.numeric(env.df$shape)
  
  # Create factor level to go with the data
  farm <- as.factor(env.df[,"FARM"])

  #### 03. Perform RDA analysis ####
  FARM_RDA <- rda(gen.imp ~ env.df$FARM)
  FARM_RDA
  
  # Perform ANOVA-like perm test for Constrained Correspondence Analysis (cca), Redundancy Analysis (rda), or...
  anova.cca(FARM_RDA, parallel=2)

  # Perform ANOVA-like perm test for rda, testing marginal effects of the terms
  anova.cca(FARM_RDA, by="margin", parallel=2)
  
  # Pull out the rsquared value
  RsquareAdj(FARM_RDA)
  
  #### 04. Plot ####
  # Set plotting variables
  bg = c("tomato1", "dodgerblue2")
  
  ### Plot RDA graph by samples ###
  RDA_plot.FN <- paste0("13_selection/rda_analysis/RDA_plot_", datatype, ".pdf")
  pdf(file = RDA_plot.FN, width = 6, height = 6)
  plot(FARM_RDA, type="n", scaling=3) # blank plot
  points(FARM_RDA, display="species", pch=20, cex=0.7, col="gray32", scaling=3)       # species = SNPs
  points(FARM_RDA, display="sites", pch=env.df$shape, cex=1, scaling=3, col=env.df$my.cols) 		# sites   = samples
  points(FARM_RDA, scaling=3, display="cn", pch = 17, col="black", cex=1.5             # cn = centroids of factor constraints
         # , pos=c(3,3,3)
         )                          # the predictors
  legend("bottomright", legend=unique(env.df$POP)
         #, bty="n"
         , col=unique(env.df$my.cols)
         , pch=16, cex=0.8
         , pt.bg=bg
         , bg = "white")
  dev.off()
  
  
  ### Identify outliers and plot RDA loadings ###
  load.rda <- summary(FARM_RDA)$species[, c("RDA1", "PC1", "PC2", "PC3")]  # extract the loadings for these axes/projections
  str(load.rda)
  
  #### Identify outliers ####
  # Make a function to identify outliers
  outliers <- function(loadings, num_sd){
    lims <- mean(loadings) + c(-1, 1) * num_sd * sd(loadings)     # find loadings +/-z sd from mean loading     
    loadings[loadings < lims[1] | loadings > lims[2]]               # locus names in these tails
  }
  
  # Set variables
  SD <- 3.5
  cand1 <- outliers(loadings = load.rda[,"RDA1"], num_sd = 3.5)
  ncand <- length(cand1)  
  
  # Plot loadings
  hist.FN <- paste0("13_selection/rda_analysis/RDA_loading_hist_", datatype, ".pdf") 
  pdf(file = hist.FN, width = 5, height = 5)
  hist_info <- hist(load.rda[, "RDA1"], main = "", xlab = "Loadings on RDA1") # Plot loadings
  text(  x = max(hist_info$mids) - max(hist_info$mids)*0.2
       , y = max(hist_info$counts) - max(hist_info$counts)*0.2
       , labels = paste0("n = ", ncand, " outliers")
       )
  abline(v = mean(load.rda[, "RDA1"]) + 3.5 * sd(load.rda[, "RDA1"]), lty = 2, col = "red")
  abline(v = mean(load.rda[, "RDA1"]) - 3.5 * sd(load.rda[, "RDA1"]), lty = 2, col = "red")
  dev.off()
  
  ## Save Results ####
  results.FN <- paste0("13_selection/rda_analysis/RDA_loading_vals_", datatype, ".csv") 
  load.rda.df <- as.data.frame(load.rda)
  load.rda.df$mname <- rownames(load.rda.df)
  write.table(x = load.rda.df, file = results.FN, append = F, quote = F, row.names = F, sep = ",")

  
  }

