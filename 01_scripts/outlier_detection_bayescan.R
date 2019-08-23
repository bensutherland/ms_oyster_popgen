# Method to detect outliers, developed by Anne-Laure Ferchaud (U Laval) and Ben Sutherland (DFO)
# Uses BayeScan methods as per http://cmpg.unibe.ch/software/BayeScan/index.html

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
#if (!require("devtools")) install.packages("devtools")
#devtools::install_github("thierrygosselin/radiator")
library(devtools)
library(radiator)

# Identify contrast files
contrasts <- list.files(path = "13_selection/sets/", pattern = ".txt")
contrasts <- contrasts[grep(pattern = "all_samples.txt", x = contrasts, invert = T)] # remove the all_samples file, is not a contrast
contrasts <- gsub(pattern = "\\.txt", replacement = "", x = contrasts) # remove the .txt

#### 1. Make Strata file for each contrast ####
# Set nulls
strata.list <- list(); samples <- NULL

# Loop
for(i in 1:length(contrasts)){
  
  # Identify samples per contrast
  contrast.FN <- paste0("13_selection/sets/", contrasts[i], ".txt")
  samples.df <- read.table(file = contrast.FN, stringsAsFactors = FALSE)
  samples.df$V2 <- gsub(samples.df$V1, pattern = "_.*", replacement = "")
  colnames(samples.df) <- c("INDIVIDUALS", "STRATA")

  # Save the strata list in the selection folder
  strata.FN <- paste0("13_selection/", contrasts[i], "_strata.txt")
  write.table(x = samples.df, file = strata.FN, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
  
}

#### 2. Prepare BayeScan input ####

##### WORKING HERE #####
for(i in 1:length(contrasts)){
  
  
  
}


genomic_converter(data = "13_selection/CHN_CHNF.recode.vcf"
                  , strata = "13_selection/CHN_CHNF_strata.txt"
                  , output = c('bayescan')
                  , filename = 'CHN_TEST_bayescan'
                  )




### EXAMPLE FROM ANNE-LAURE
#### FARM SELECTION ####
FRA_FRAF = genomic_converter(data = "13_selection/FRA_FRAF.recode.vcf",
                             strata = "13_selection/FRA_FRAF_strata.txt",
                             output =c('bayescan'),
                             monomorphic.out = FALSE,
                             common.markers = FALSE,
                             vcf.metadata = TRUE,
                             filename = 'FRA_FRAF_bayescan'
)
#### END EXAMPLE ####



# Identify files
files <- list.files(path = "13_selection/", pattern = "recode.vcf", full.names = T)





####  STEP 1: make input and run BAYESCAN ####


# e.g. ##
# set nulls
pcadapt.obj <- NULL; pcadapt.FN <- NULL; pcadapt.out <- NULL
pcadapt.out.list <- list(); obj.name <- NULL; pcadapt.obj.list <- list()

for(i in 1:length(files)){
  
  # Select file
  pcadapt.FN <- files[i]
  print(paste0("Working on ", pcadapt.FN))
  
  # Create object name
  obj.name <- gsub(pattern = ".*//", replacement = "", x = files[i])
  obj.name <- gsub(pattern = "\\..*", replacement = "", x = obj.name)
  
  # Read in data
  pcadapt.obj <- read.pcadapt(input = pcadapt.FN, type = "vcf")
  
  # Save input
  pcadapt.obj.list[[obj.name]] <- pcadapt.obj
  
  # Run pcadapt
  print("Running pcadapt with K = 20")
  pcadapt.out <- pcadapt(input = pcadapt.obj, K = 20)
  
  # Save output
  pcadapt.out.list[[obj.name]] <- pcadapt.out
  
}
# end e.g. ##

# Use genomic_converter to prep Bayescan input per contrast based on the recoded VCF and the strata file (to be built)






#run in terminal with pr_odd = 100
./BayeScan2.1_macos64bits ../../FARM_SELECTION/FRA_FRAF_bayescan.txt ../../FARM_SELECTION/FRA_FRAF_bayescan -threads 5 -pr_odds 100


CHN_CHNF = genomic_converter(data = "CHN_CHNF.recode.vcf",
                             strata = "CHN_CHNF_strata.txt",
                             output =c('bayescan'),
                             monomorphic.out = FALSE,
                             common.markers = FALSE,
                             vcf.metadata = TRUE,
                             filename = 'CHN_CHNF_bayescan'
)
#run in terminalwith pr_odd = 100
./BayeScan2.1_macos64bits ../../FARM_SELECTION/CHN_CHNF_bayescan.txt ../../FARM_SELECTION/CHN_CHNF_bayescan -threads 5 -pr_odds 100


PEN_PENF = genomic_converter(data = "PEN_PENF.recode.vcf",
                             strata = "PEN_PENF_strata.txt",
                             output =c('bayescan'),
                             monomorphic.out = FALSE,
                             common.markers = FALSE,
                             vcf.metadata = TRUE,
                             filename = 'PEN_PENF_bayescan'
)
#run in terminal with pr_odd = 100
./BayeScan2.1_macos64bits ../../FARM_SELECTION/PEN_PENF_bayescan.txt ../../FARM_SELECTION/PEN_PENF_bayescan -threads 5 -pr_odds 100



#### DOMESTIC SELECTION ####

setwd('/Volumes/drobo2/Anne-Laure/PACIFIC_OYSTER/DOMESTIC_SELECTION')
DBP_PEN = genomic_converter(data = "DBP_PEN.recode.vcf",
                            strata = "DBP_PEN_strata.txt",
                            output =c('bayescan'),
                            monomorphic.out = FALSE,
                            common.markers = FALSE,
                            vcf.metadata = TRUE,
                            filename = 'DBP_PEN_bayescan'
)
#run in terminal with pr_odd = 100
./BayeScan2.1_macos64bits ../../DOMESTIC_SELECTION/DBP_PEN_bayescan.txt -o ../../DOMESTIC_SELECTION/DBP_PEN_bayescan -threads 5 -pr_odds 100

ROS_PIP = genomic_converter(data = "ROS_PIP.recode.vcf",
                            strata = "ROS_PIP_strata.txt",
                            output =c('bayescan'),
                            monomorphic.out = FALSE,
                            common.markers = FALSE,
                            vcf.metadata = TRUE,
                            filename = 'ROS_PIP_bayescan'
)
#run in terminal with pr_odd = 100
./BayeScan2.1_macos64bits ../../DOMESTIC_SELECTION/ROS_PIP_bayescan.txt -o ../../DOMESTIC_SELECTION/ROS_PIP_bayescan -threads 5 -pr_odds 100

GUR_FRA = genomic_converter(data = "GUR_FRA.recode.vcf",
                            strata = "GUR_FRA_strata.txt",
                            output =c('bayescan'),
                            monomorphic.out = FALSE,
                            common.markers = FALSE,
                            vcf.metadata = TRUE,
                            filename = 'GUR_FRA_bayescan'
)
#run in terminal with pr_odd = 100
./BayeScan2.1_macos64bits ../../DOMESTIC_SELECTION/GUR_FRA_bayescan.txt -o ../../DOMESTIC_SELECTION/GUR_FRA_bayescan -threads 5 -pr_odds 100

QDC_CHN = genomic_converter(data = "QDC_CHN.recode.vcf",
                            strata = "QDC_CHN_strata.txt",
                            output =c('bayescan'),
                            monomorphic.out = FALSE,
                            common.markers = FALSE,
                            vcf.metadata = TRUE,
                            filename = 'QDC_CHN_bayescan'
)

#run in terminal with pr_odd = 100
./BayeScan2.1_macos64bits ../../DOMESTIC_SELECTION/QDC_CHN_bayescan.txt -o ../../DOMESTIC_SELECTION/QDC_CHN_bayescan -threads 5 -pr_odds 100

RSC_CHN = genomic_converter(data = "RSC_CHN.recode.vcf",
                            strata = "RSC_CHN_strata.txt",
                            output =c('bayescan'),
                            monomorphic.out = FALSE,
                            common.markers = FALSE,
                            vcf.metadata = TRUE,
                            filename = 'RSC_CHN_bayescan'
)

#run Bayescan in Terminal with pr_odd = 100
./BayeScan2.1_macos64bits ../../DOMESTIC_SELECTION/RSC_CHN_bayescan.txt -o ../../DOMESTIC_SELECTION/RSC_CHN_bayescan -threads 5 -pr_odds 100





####  STEP 2:  work with BAYESCAN results ####

source ("/Users/ALF/Documents/10-PROGRAMS/BayeScan2.1/R functions/plot_R.r")

#### FARM SELECTION ###
setwd('/Users/ALF/Dropbox/03-QUEBEC/16-PACIFIC_OYSTER/BAYESCAN/FARM_SELECTION')

#example w/ CHN / CHNF

#CHN_CHNF results####
bayescan_plot_CHN_CHNF = plot_bayescan("CHN_CHNF_bayescan_fst.txt", FDR = 0.01)

bayescan_res_CHN_CHNF= read.table("CHN_CHNF_bayescan_fst.txt")
SNPb_CHN_CHNF=read.table("CHN_CHNF_bayescan_markers_dictionary.tsv",header=T)
bayescan_CHN_CHNF=cbind(SNPb_CHN_CHNF$MARKERS, bayescan_res_CHN_CHNF)
colnames(bayescan_CHN_CHNF)=c("SNP","POST_PROB","LOG10_PO","Q_VALUE","ALPHA","FST")

# POST_PROB = 1 & Q_VALUE = 0 == 0.0001 
attach(bayescan_CHN_CHNF)
class(bayescan_CHN_CHNF$Q_VALUE)
bayescan_CHN_CHNF$Q_VALUE <- as.numeric(bayescan_CHN_CHNF$Q_VALUE)
bayescan_CHN_CHNF[bayescan_CHN_CHNF$Q_VALUE<=0.0001,"Q_VALUE"]=0.0001

bayescan_CHN_CHNF$POST_PROB <- (round(bayescan_CHN_CHNF$POST_PROB, 4))
bayescan_CHN_CHNF$LOG10_PO <- (round(bayescan_CHN_CHNF$LOG10_PO, 4))
bayescan_CHN_CHNF$Q_VALUE <- (round(bayescan_CHN_CHNF$Q_VALUE, 4))
bayescan_CHN_CHNF$ALPHA <- (round(bayescan_CHN_CHNF$ALPHA, 4))
bayescan_CHN_CHNF$FST <- (round(bayescan_CHN_CHNF$FST, 6))

# ADD a column for selection grouping
bayescan_CHN_CHNF$SELECTION <- ifelse(bayescan_CHN_CHNF$ALPHA>=0&bayescan_CHN_CHNF$Q_VALUE<=0.05,"diversifying",ifelse(bayescan_CHN_CHNF$ALPHA>=0&bayescan_CHN_CHNF$Q_VALUE>0.05,"neutral","balancing"))
bayescan_CHN_CHNF$SELECTION<- factor(bayescan_CHN_CHNF$SELECTION)
levels(bayescan_CHN_CHNF$SELECTION)

positive_CHN_CHNF <- bayescan_CHN_CHNF[bayescan_CHN_CHNF$SELECTION=="diversifying",]
neutral_CHN_CHNF <- bayescan_CHN_CHNF[bayescan_CHN_CHNF$SELECTION=="neutral",]
balancing_CHN_CHNF <- bayescan_CHN_CHNF[bayescan_CHN_CHNF$SELECTION=="balancing",]
xtabs(data=bayescan_CHN_CHNF, ~SELECTION)
write.table(neutral_CHN_CHNF, "neutral_CHN_CHNF.txt", row.names=F, quote=F)
write.table(balancing_CHN_CHNF, "balancing_CHN_CHNF.txt", row.names=F, quote=F)
write.table(positive_CHN_CHNF, "positive_CHN_CHNF", row.names=F, quote=F)
detach(bayescan_CHN_CHNF)
