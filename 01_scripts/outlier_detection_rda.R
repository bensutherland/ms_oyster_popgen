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

for(i in 1:length(datatypes)){
  
  # Select your datatype
  datatype <- datatypes[i]
  
  ##  Read in genotype data ##
  geno.FN <- paste0("13_selection/rda_analysis/", datatype, "_geno.012")
  geno <- read.table(file = geno.FN)
  ind <- geno[,1]
  gen <- as.matrix(geno[,2:ncol(geno)])
  
  # Identify NAs
  sum(is.na(gen)) # no NAs
  gen[gen == -1] <- NA
  sum(is.na(gen)) # NAs now in data
  
  # Impute
  gen.imp <- apply(X = gen, MARGIN = 2, function(x)replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
  # work over columns (replace the missing data with the most common genotype in the dataset (?))
  sum(is.na(gen.imp))
  
  #### Make env data ####
  # e.g. in the following format:
  # INDIVIDUALS	POP	FARM
  # CHN_1551	CHN	WILD
  # CHN_1552	CHN	WILD
  # CHN_1553	CHN	WILD
  # CHN_1554	CHN	WILD
  
  # Read in the samples in this comparison
  env.FN <- paste0("13_selection/rda_analysis/", datatype, ".txt")
  env <- read.table(env.FN)
  env.df <- as.data.frame(env)
  env.df$V2 <- gsub(pattern = "_.*", replacement = "", x =  env.df$V1)
  head(env.df)
  
  if(datatype=="FARM"){
    # Identify the farm pops
    selection_pops <- c("CHNF", "PENF", "FRAF")
  }else if(datatype=="DOMESTICATION"){
    # Identify the farm pops
    selection_pops <- c("DPB", "GUR", "QDC", "ROS", "RSC")
  }
  
  # Identify which are the selection pops
  env.df$V3[env.df$V2 %in% selection_pops] <- "FARM"
  env.df$V3[!env.df$V2 %in% selection_pops] <- "WILD" # find opposite of the above
  
  colnames(env.df) <- c("INDIVIDUALS", "POP", "FARM")
  
  # OLD METHOD:
  # env.df$V3[grep(pattern = "F_", x = env.df$V1)] <- "FARM"
  # env.df$V3[grep(pattern = "F_", x = env.df$V1, invert = T)] <- "WILD"
  # END OLD METHOD
  
  # Inspect
  # env.df
  head(env.df)
  tail(env.df)
  
  # Create factor level to go with the data
  farm <- as.factor(env.df[,"FARM"])

  
  #### RDA Analysis ####
  FARM_RDA <- rda(gen.imp ~ farm)
  FARM_RDA
  
  # Perform ANOVA-like perm test for Constrained Correspondence Analysis (cca), Redundancy Analysis (rda), or...
  anova.cca(FARM_RDA, parallel=2)

  # Perform ANOVA-like perm test for rda, testing marginal effects of the terms
  anova.cca(FARM_RDA, by="margin", parallel=2)
  
  # 
  RsquareAdj(FARM_RDA)
  
  
  #### Plotting ####
  bg = c("tomato1", "dodgerblue2")
  # [1] "#7FC97F" "#BEAED4" "#FDC086"
  
  ### RDA graph by pop ###
  # pdf("POP_RDA1-2_FARM.pdf", 9, 9)
  plot(FARM_RDA, type="n", scaling=3) # blank plot
  
  # add points
  points(FARM_RDA, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
  points(FARM_RDA, display="sites", pch=17, cex=1, scaling=3, col=bg[farm]) 		# the fish
  points(FARM_RDA, scaling=3, display="cn", pch =17, col="black", cex=1.5
         # , pos=c(3,3,3)
         )                          # the predictors
  legend("bottomright", inset=c(0, 0), legend=levels(farm), bty="n", col=bg, pch=17, cex=1, pt.bg=bg)
  # dev.off()
  
  ### Question: 
  ## What does the PCA vs RDA mean? Need more detail on interp. 
  
  ### RDA graph (farm status: symbol; pop: colour)
  #Create colors
  library(RColorBrewer)
  n <- length(unique(x = env.df$POP)) # the number unique colours needed
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  bg = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  bg = bg[1:n]
  
  # bc = c("17","18") # not sure what this is.. (bjgs)
  
  pop <- as.factor(env.df[,2])
  
  ### RDA graph by pop with colours
  #pdf("POP_RDA1-2_FARM_by_pop.pdf", 9, 9)
  plot(FARM_RDA, type="n", scaling=3)
  
  # add points
  points(FARM_RDA, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
  points(FARM_RDA, display="sites", pch=c(17, 18)[as.factor(farm)], cex=1, scaling=3, col=bg[pop]) 		# the fish
  legend("bottomleft", 
         legend = unique(pop), 
         col = bg, 
         pch = 17,  # TO FIX DYNAMIC *#TODO* BJGS
         bty = "n", 
         pt.cex = 2, 
         cex = 1.2, 
         text.col = "black", 
         horiz = F , 
         inset = c(0.01, 0.01))
  # dev.off()
  
  
  ##### HERE TODAY #####
  
  # RDA GRAPH by SNPs
  load.rda <- summary(FARM_RDA)$species[,1:4] 
  hist(load.rda[,1], main="Loadings on RDA1")
  
  outliers <- function(x,z){
    lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
    x[x < lims[1] | x > lims[2]]               # locus names in these tails
  }
  
  #Choose SD
  SD=3.5
  cand1 <- outliers(load.rda[,1], 3.5) #46
  ncand <- length(cand1)  #46
  
  cand <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
  colnames(cand) <- c("axis","snp","loading")
  # write.table(cand, file="listeCandidat_with_farm.txt", row.names = F)
  # OUTPUT
  
}

